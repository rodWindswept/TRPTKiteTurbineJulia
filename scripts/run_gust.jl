# scripts/run_gust.jl
# Interactive dashboard with a raised-cosine (Hann-window) gust event.
# Baseline 11 m/s with a gust to 17 m/s peaking at t=60 s (t=40–80 s window).
# Demonstrates TRPT response to a transient over-speed event.
#
# Usage:
#   julia --project=. scripts/run_gust.jl

using DifferentialEquations, GLMakie, Statistics

include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/dynamics.jl")
include("../src/geometry.jl")
include("../src/force_analysis.jl")
include("../src/visualization.jl")

mkpath("output")

const SIM_DURATION = 120.0

p = params_10kw()

# Hann-window gust: baseline 11 m/s → peak 17 m/s over t=40–80 s
wind_fn = gust_event(11.0, 17.0, 40.0, 80.0)

println("Gust event: 11 m/s baseline, 17 m/s peak at t=60 s (window t=40–80 s)")
println("Solving $(SIM_DURATION) s ODE...")

sol = solve(ODEProblem(trpt_ode_wind!, [0.0, 1.0], (0.0, SIM_DURATION), (p, wind_fn)),
            Tsit5(); reltol=1e-6, abstol=1e-6, saveat=1/30)
println("  $(length(sol.t)) frames  retcode: $(sol.retcode)")

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
println("Max V_hub: $(round(maximum(v_hub_ts), digits=2)) m/s  at t=$(sol.t[argmax(v_hub_ts)]) s")

fig, _ = build_trpt_scene(p, traj)
display(fig)
readline()
