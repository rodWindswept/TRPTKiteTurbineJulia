# scripts/run_wind_ramp.jl
# Interactive dashboard with a linear wind ramp: 6 → 14 m/s over t = 20–100 s.
# Demonstrates rotor response to a gradual increase in wind speed.
#
# Usage:
#   julia --project=. scripts/run_wind_ramp.jl

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

p = params_10kw()

# Wind function: reference speed ramps from 6 → 14 m/s between t=20 and t=100 s
wind_fn = wind_ramp(6.0, 14.0, 20.0, 100.0)

println("Wind ramp: 6 → 14 m/s over t=20–100 s")
println("Solving $(SIM_DURATION) s ODE...")

sol = solve(ODEProblem(trpt_ode_wind!, [0.0, 1.0], (0.0, SIM_DURATION), (p, wind_fn)),
            Tsit5(); reltol=1e-6, abstol=1e-6, saveat=1/30)
println("  $(length(sol.t)) frames  retcode: $(sol.retcode)")

# Build v_hub time series at hub altitude for each frame
h_hub    = hub_altitude(p.tether_length, p.elevation_angle)
v_hub_ts = [wind_at_altitude(wind_fn(t), p.h_ref, h_hub) for t in sol.t]

traj = (
    t         = sol.t,
    alpha_tot = sol[1, :],
    omega     = sol[2, :],
    power_kw  = [instantaneous_power(p, [sol[1,i], sol[2,i]]) / 1000.0
                 for i in eachindex(sol.t)],
    v_hub     = v_hub_ts,
)

println("Mean $(round(mean(traj.power_kw), digits=2)) kW  |  Peak $(round(maximum(traj.power_kw), digits=2)) kW")

fig, _ = build_trpt_scene(p, traj)
display(fig)
readline()
