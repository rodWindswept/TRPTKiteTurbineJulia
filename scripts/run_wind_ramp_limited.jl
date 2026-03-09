# scripts/run_wind_ramp_limited.jl
# Wind ramp 6 → 23 m/s over 120 s with elevation-angle power limiter active.
# Demonstrates MPPT region (β flat at 23°) transitioning to limited region
# (β rising to shed excess power) as wind passes rated speed (~11 m/s at ~55 s).
# Compare with run_wind_ramp.jl which uses the 2-state MPPT ODE (no limiter).
#
# Usage:
#   julia --project=. scripts/run_wind_ramp_limited.jl

using DifferentialEquations, GLMakie, Statistics

include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/aerodynamics.jl")
include("../src/dynamics.jl")
include("../src/geometry.jl")
include("../src/force_analysis.jl")
include("../src/visualization.jl")

mkpath("output")

const SIM_DURATION = 120.0
const V_START      = 6.0     # m/s
const V_END        = 23.0    # m/s

p       = params_10kw()
wind_fn = wind_ramp(V_START, V_END, 20.0, 100.0)

# Warm-start: optimal TSR at the starting wind speed
h_hub  = hub_altitude(p.tether_length, p.β_min)
v_hub0 = wind_at_altitude(V_START, p.h_ref, h_hub)
ω0     = 4.1 * v_hub0 / p.rotor_radius
u0     = [0.0, ω0, p.β_min]

println("Wind ramp: $(V_START) → $(V_END) m/s over t=20–100 s  (elevation limiter active)")
println("Limiter: P_rated=$(p.p_rated_w/1000) kW, β_min=$(round(rad2deg(p.β_min),digits=1))°, β_max=$(round(rad2deg(p.β_max),digits=1))°")
println("Solving $(SIM_DURATION) s ODE...")

sol = solve(ODEProblem(trpt_ode_wind_limited!, u0, (0.0, SIM_DURATION), (p, wind_fn)),
            Tsit5(); reltol=1e-6, abstol=1e-6, saveat=1/30)
println("  $(length(sol.t)) frames  retcode: $(sol.retcode)")

v_hub_ts = [wind_at_altitude(wind_fn(t), p.h_ref,
             hub_altitude(p.tether_length, sol[3, i]))
            for (i, t) in enumerate(sol.t)]

traj = (
    t         = sol.t,
    alpha_tot = sol[1, :],
    omega     = sol[2, :],
    power_kw  = [instantaneous_power(p, [sol[1,i], sol[2,i]], sol[3,i]) / 1000.0
                 for i in eachindex(sol.t)],
    v_hub     = v_hub_ts,
    beta      = rad2deg.(sol[3, :]),
)

println("Peak power: $(round(maximum(traj.power_kw), digits=2)) kW")
println("β range:    $(round(minimum(traj.beta), digits=1))° – $(round(maximum(traj.beta), digits=1))°")
println("Mean power (last 30 s): $(round(mean(traj.power_kw[end-round(Int,30*30):end]), digits=2)) kW")

fig, _ = build_trpt_scene(p, traj)
display(fig)
readline()
