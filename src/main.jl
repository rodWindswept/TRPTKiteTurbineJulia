# src/main.jl
# Simulation orchestration for the TRPT Kite Turbine Simulator.
# Runs a 120 s simulation, exports telemetry CSV, and launches GLMakie animation.
#
# Usage:
#   julia --project=. src/main.jl           # default 10 kW
#   julia --project=. src/main.jl 50        # 50 kW

using DifferentialEquations
using DataFrames
using CSV
using GLMakie
using Statistics

include("parameters.jl")
include("wind_profile.jl")
include("aerodynamics.jl")
include("dynamics.jl")
include("geometry.jl")
include("force_analysis.jl")
include("visualization.jl")

const RATED_KW       = length(ARGS) > 0 ? parse(Float64, ARGS[1]) : 10.0
const SIM_DURATION_S = 120.0

function select_params(rated_kw::Float64)::SystemParams
    if rated_kw ≈ 10.0
        return params_10kw()
    elseif rated_kw ≈ 50.0
        return params_50kw()
    else
        println("Custom rating $rated_kw kW — scaling from 10 kW baseline.")
        return mass_scale(params_10kw(), 10.0, rated_kw)
    end
end

function main()
    mkpath("output")

    println("="^60)
    println("TRPT Kite Turbine Simulator — $RATED_KW kW configuration")
    println("="^60)

    p     = select_params(RATED_KW)
    h_hub = hub_altitude(p.tether_length, p.elevation_angle)
    v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h_hub)
    n_seg = p.n_rings + 1

    # Total inertia (same formula as trpt_ode!)
    r_top_i = p.trpt_hub_radius
    r_bot_i = 2.0 * p.tether_length * p.trpt_rL_ratio / n_seg - r_top_i
    I_rings = sum(p.m_ring * (r_bot_i + i / p.n_rings * (r_top_i - r_bot_i))^2
                  for i in 1:p.n_rings)
    I_total = p.n_blades * p.m_blade * p.rotor_radius^2 + I_rings + p.i_pto

    println("Hub altitude   : $(round(h_hub, digits=1)) m")
    println("Wind at hub    : $(round(v_hub, digits=2)) m/s  (ref $(p.v_wind_ref) m/s @ $(p.h_ref) m)")
    println("Rotor radius   : $(p.rotor_radius) m")
    println("n_lines / rings: $(p.n_lines) / $(p.n_rings)")
    println("Total inertia  : $(round(I_total, digits=2)) kg·m²")
    println()

    # ── ODE solve ─────────────────────────────────────────────────────────────
    println("Solving $(SIM_DURATION_S) s ODE (Tsit5)...")
    u0   = [0.0, 1.0]
    prob = ODEProblem(trpt_ode!, u0, (0.0, SIM_DURATION_S), p)
    sol  = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=1/30)

    if sol.retcode != ReturnCode.Success
        error("ODE solver failed: retcode=$(sol.retcode)")
    end
    println("Solved $(length(sol.t)) time points. retcode: $(sol.retcode)")

    # ── Telemetry arrays ──────────────────────────────────────────────────────
    power_kw = [instantaneous_power(p, [sol[1, i], sol[2, i]]) / 1000.0
                for i in eachindex(sol.t)]

    traj = (
        t         = sol.t,
        alpha_tot = sol[1, :],
        omega     = sol[2, :],
        power_kw  = power_kw,
        v_hub     = v_hub,
    )

    # ── Force scan over full trajectory ───────────────────────────────────────
    u_frames  = [[sol[1, i], sol[2, i]] for i in eachindex(sol.t)]
    v_hub_all = fill(v_hub, length(sol.t))
    T_max_run, C_max_run = run_force_scan(p, u_frames, v_hub_all)
    force_states = [element_forces(p, u_frames[i], v_hub_all[i]) for i in eachindex(sol.t)]

    # ── CSV export ────────────────────────────────────────────────────────────
    df = DataFrame(
        time_s            = traj.t,
        total_twist_deg   = rad2deg.(traj.alpha_tot),
        segment_twist_deg = rad2deg.(traj.alpha_tot ./ n_seg),
        rotor_speed_rads  = traj.omega,
        hub_wind_speed_ms = fill(v_hub, length(traj.t)),
        power_kw          = traj.power_kw,
        altitude_m        = fill(h_hub, length(traj.t)),
        tau_aero          = [fs.tau_aero        for fs in force_states],
        tau_drag          = [fs.tau_drag        for fs in force_states],
        tau_transmitted   = [fs.tau_transmitted for fs in force_states],
        T_max_run         = fill(T_max_run, length(sol.t)),
        C_max_run         = fill(C_max_run, length(sol.t)),
    )

    for seg in 1:n_seg
        df[!, Symbol("T_seg_$(seg)")] = [fs.tether_tension[seg] for fs in force_states]
    end
    for ring in 1:p.n_rings
        df[!, Symbol("C_ring_$(ring)")] = [fs.ring_compression[ring] for fs in force_states]
    end

    csv_path = "output/telemetry_summary.csv"
    CSV.write(csv_path, df)
    println("Telemetry exported → $csv_path  ($(nrow(df)) rows)")

    # ── Summary statistics ────────────────────────────────────────────────────
    peak_power    = maximum(df.power_kw)
    mean_power    = mean(df.power_kw)
    max_twist_deg = maximum(abs.(df.segment_twist_deg))
    n_collapse    = sum(abs.(df.segment_twist_deg) .>= 0.95 * 180.0)

    println()
    println("── Flight summary ──────────────────────────────")
    println("Peak ground power   : $(round(peak_power, digits=2)) kW")
    println("Mean ground power   : $(round(mean_power, digits=2)) kW")
    println("Max segment twist   : $(round(max_twist_deg, digits=1))°  (limit: 171°)")
    println("Collapse events     : $n_collapse")
    println("────────────────────────────────────────────────")

    # ── GLMakie animation ─────────────────────────────────────────────────────
    println("\nLaunching 3D visualization...")
    fig, _ = build_trpt_scene(p, traj)
    save("output/main_simulation_frame.png", fig)
    println("Saved initial frame → output/main_simulation_frame.png")
    display(fig)
    println("Press Enter to exit after viewing the visualization.")
    readline()
end

main()
