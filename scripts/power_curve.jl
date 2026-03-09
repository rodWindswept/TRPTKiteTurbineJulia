# scripts/power_curve.jl
# Sweep wind speeds from 4 to 16 m/s (1 m/s steps) and compute steady-state
# ground power for both 10 kW and 50 kW configurations.
# Outputs a printed table, output/power_curve.csv, and output/power_curve.png.
#
# Usage:
#   julia --project=. scripts/power_curve.jl

using DifferentialEquations, DataFrames, CSV, Statistics, Plots

include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/aerodynamics.jl")
include("../src/dynamics.jl")

mkpath("output")

const WIND_SPEEDS   = 4.0:1.0:23.0   # m/s at reference height p.h_ref (extended to 23 m/s)
const SETTLE_TIME   = 120.0           # allow 120 s for steady state (limiter needs longer to settle)
const SAMPLE_START  = 100.0           # average power over final 20 s

function steady_power_kw(p::SystemParams, v_ref::Float64)::Float64
    # Warm-start near optimal TSR so MPPT controller locks on quickly
    h_hub = hub_altitude(p.tether_length, p.elevation_angle)
    v_hub = wind_at_altitude(v_ref, p.h_ref, h_hub)
    ω0 = 4.1 * v_hub / p.rotor_radius   # λ_opt × v_hub / R
    wind_fn = steady_wind(v_ref)
    sol = solve(ODEProblem(trpt_ode_wind!, [0.0, ω0], (0.0, SETTLE_TIME), (p, wind_fn)),
                Tsit5(); reltol=1e-6, abstol=1e-6, saveat=0.1)
    sol.retcode == ReturnCode.Success || return NaN
    # Average power during final 20 s (steady state)
    sample_idx = findall(t -> t >= SAMPLE_START, sol.t)
    isempty(sample_idx) && return NaN
    mean_power = mean(instantaneous_power(p, [sol[1,i], sol[2,i]])
                      for i in sample_idx)
    return mean_power / 1000.0
end

function steady_power_kw_limited(p::SystemParams, v_ref::Float64)::Float64
    # 3-state limiter ODE: β starts at β_min and rises to shed excess power above rated
    h_hub = hub_altitude(p.tether_length, p.β_min)
    v_hub = wind_at_altitude(v_ref, p.h_ref, h_hub)
    ω0 = 4.1 * v_hub / p.rotor_radius
    wind_fn = steady_wind(v_ref)
    sol = solve(ODEProblem(trpt_ode_wind_limited!, [0.0, ω0, p.β_min], (0.0, SETTLE_TIME), (p, wind_fn)),
                Tsit5(); reltol=1e-6, abstol=1e-6, saveat=0.1)
    sol.retcode == ReturnCode.Success || return NaN
    sample_idx = findall(t -> t >= SAMPLE_START, sol.t)
    isempty(sample_idx) && return NaN
    mean_power = mean(instantaneous_power(p, [sol[1,i], sol[2,i]], sol[3,i])
                      for i in sample_idx)
    return mean_power / 1000.0
end

configs = [("10 kW", params_10kw()), ("50 kW", params_50kw())]

println("="^60)
println("TRPT Power Curve — steady-state output vs wind speed")
println("="^60)

all_rows = DataFrame(v_ref_ms = Float64[], config = String[], mode = String[], power_kw = Float64[])

for (label, p) in configs
    println("\n$label configuration (rotor $(p.rotor_radius) m):")
    println("  V_ref (m/s) │ MPPT (kW) │ Limited (kW)")
    println("  ────────────┼───────────┼─────────────")
    for v in WIND_SPEEDS
        pw_mppt = steady_power_kw(p, Float64(v))
        pw_lim  = steady_power_kw_limited(p, Float64(v))
        println("  $(lpad(v, 10)) │ $(rpad(round(pw_mppt, digits=2), 9)) │ $(round(pw_lim, digits=2))")
        push!(all_rows, (v_ref_ms=Float64(v), config=label, mode="MPPT",    power_kw=pw_mppt))
        push!(all_rows, (v_ref_ms=Float64(v), config=label, mode="Limited", power_kw=pw_lim))
    end
end

csv_path = "output/power_curve.csv"
CSV.write(csv_path, all_rows)
println("\nSaved → $csv_path  ($(nrow(all_rows)) rows)")

# Plot all 4 curves on one chart
plt = plot(;
    xlabel = "Wind speed at hub (m/s)",
    ylabel = "Ground power (kW)",
    title  = "TRPT Kite Turbine — Power Curve  (MPPT vs Elevation Limiter)",
    legend = :topleft,
    grid   = true,
    size   = (900, 550),
)
styles = Dict("MPPT" => :solid, "Limited" => :dash)
colors = Dict("10 kW" => :steelblue, "50 kW" => :firebrick)
for (label, p) in configs
    for mode in ["MPPT", "Limited"]
        rows = filter(r -> r.config == label && r.mode == mode, all_rows)
        plot!(plt, rows.v_ref_ms, rows.power_kw;
              label = "$label $mode  (R=$(round(p.rotor_radius, digits=1)) m)",
              lw = 2, marker = :circle, markersize = 4,
              linestyle = styles[mode], color = colors[label])
    end
end

png_path = "output/power_curve.png"
savefig(plt, png_path)
println("Saved → $png_path")
