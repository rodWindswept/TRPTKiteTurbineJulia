# scripts/viz_smoke_test.jl — GLMakie scene smoke test
using GLMakie, DifferentialEquations
include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/dynamics.jl")
include("../src/visualization.jl")

p  = params_10kw()
u0 = [0.0, 1.0]
prob = ODEProblem(trpt_ode!, u0, (0.0, 30.0), p)
sol  = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=1/30)

n_seg = p.n_rings + 1
h_hub = hub_altitude(p.tether_length, p.elevation_angle)
v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h_hub)

traj = (
    t         = sol.t,
    alpha_tot = sol[1, :],
    omega     = sol[2, :],
    power_kw  = [instantaneous_power(p, [sol[1, i], sol[2, i]]) / 1000.0
                 for i in 1:length(sol.t)],
    v_hub     = v_hub,
)

fig, _ = build_trpt_scene(p, traj)
save("output/viz_smoke_test.png", fig)
println("Smoke test passed — saved output/viz_smoke_test.png")
