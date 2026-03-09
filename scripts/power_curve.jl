# scripts/power_curve.jl
# Sweep wind speeds from 4 to 16 m/s (1 m/s steps) and compute steady-state
# ground power for both 10 kW and 50 kW configurations.
# Outputs a printed table and saves output/power_curve.csv.
#
# Usage:
#   julia --project=. scripts/power_curve.jl

using DifferentialEquations, DataFrames, CSV, Statistics

include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/aerodynamics.jl")
include("../src/dynamics.jl")

mkpath("output")

const WIND_SPEEDS   = 4.0:1.0:16.0   # m/s at reference height p.h_ref
const SETTLE_TIME   = 60.0            # allow 60 s for steady state
const SAMPLE_START  = 50.0            # average power over final 10 s

function steady_power_kw(p::SystemParams, v_ref::Float64)::Float64
    # Warm-start near optimal TSR so MPPT controller locks on quickly
    h_hub = hub_altitude(p.tether_length, p.elevation_angle)
    v_hub = wind_at_altitude(v_ref, p.h_ref, h_hub)
    ω0 = 4.1 * v_hub / p.rotor_radius   # λ_opt × v_hub / R
    wind_fn = steady_wind(v_ref)
    sol = solve(ODEProblem(trpt_ode_wind!, [0.0, ω0], (0.0, SETTLE_TIME), (p, wind_fn)),
                Tsit5(); reltol=1e-6, abstol=1e-6, saveat=0.1)
    sol.retcode == ReturnCode.Success || return NaN
    # Average power during final 10 s (steady state)
    sample_idx = findall(t -> t >= SAMPLE_START, sol.t)
    isempty(sample_idx) && return NaN
    mean_power = mean(instantaneous_power(p, [sol[1,i], sol[2,i]])
                      for i in sample_idx)
    return mean_power / 1000.0
end

configs = [("10 kW", params_10kw()), ("50 kW", params_50kw())]

println("="^60)
println("TRPT Power Curve — steady-state output vs wind speed")
println("="^60)

all_rows = DataFrame(v_ref_ms = Float64[], config = String[], power_kw = Float64[])

for (label, p) in configs
    println("\n$label configuration (rotor $(p.rotor_radius) m):")
    println("  V_ref (m/s) │ Power (kW)")
    println("  ────────────┼───────────")
    for v in WIND_SPEEDS
        pw = steady_power_kw(p, Float64(v))
        println("  $(lpad(v, 10)) │ $(rpad(round(pw, digits=2), 10))")
        push!(all_rows, (v_ref_ms=Float64(v), config=label, power_kw=pw))
    end
end

csv_path = "output/power_curve.csv"
CSV.write(csv_path, all_rows)
println("\nSaved → $csv_path  ($(nrow(all_rows)) rows)")
