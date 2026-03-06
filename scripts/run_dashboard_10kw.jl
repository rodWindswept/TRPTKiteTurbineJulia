# scripts/run_dashboard_10kw.jl
# Interactive 3D GLMakie dashboard — 10 kW prototype at rated wind.
#
# Usage:
#   julia --project=. scripts/run_dashboard_10kw.jl

using DifferentialEquations, DataFrames, CSV, GLMakie, Statistics

include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/aerodynamics.jl")
include("../src/dynamics.jl")
include("../src/geometry.jl")
include("../src/force_analysis.jl")
include("../src/visualization.jl")

mkpath("output")

p     = params_10kw()
h_hub = hub_altitude(p.tether_length, p.elevation_angle)
v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h_hub)
n_seg = p.n_rings + 1

println("10 kW config — hub $(round(h_hub, digits=1)) m, V_hub $(round(v_hub, digits=2)) m/s")
println("Solving 120 s ODE...")

sol = solve(ODEProblem(trpt_ode!, [0.0, 1.0], (0.0, 120.0), p),
            Tsit5(); reltol=1e-6, abstol=1e-6, saveat=1/30)

traj = (
    t         = sol.t,
    alpha_tot = sol[1, :],
    omega     = sol[2, :],
    power_kw  = [instantaneous_power(p, [sol[1,i], sol[2,i]]) / 1000.0
                 for i in eachindex(sol.t)],
    v_hub     = v_hub,
)

println("Mean $(round(mean(traj.power_kw), digits=2)) kW  |  Peak $(round(maximum(traj.power_kw), digits=2)) kW")

fig, _ = build_trpt_scene(p, traj)
display(fig)
readline()
