# scripts/run_turbulent.jl
# Interactive dashboard with IEC 61400-1 Class A turbulence (TI = 15 %).
# Uses a first-order Markov (AR(1)) wind model with L = 340 m integral length
# scale. The wind time series is pre-computed then linearly interpolated.
#
# Usage:
#   julia --project=. scripts/run_turbulent.jl

using DifferentialEquations, GLMakie, Statistics

include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/dynamics.jl")
include("../src/visualization.jl")

mkpath("output")

const SIM_DURATION  = 120.0
const V_MEAN        = 11.0    # m/s (rated wind speed)
const TURBULENCE_I  = 0.15    # 15 % TI — IEC Class A onshore
const RNG_SEED      = 42

p       = params_10kw()
wind_fn = turbulent_wind(V_MEAN, TURBULENCE_I, SIM_DURATION; rng_seed=RNG_SEED)

println("Turbulent wind: mean=$(V_MEAN) m/s, TI=$(round(TURBULENCE_I*100))%, seed=$(RNG_SEED)")
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

fig, _ = build_trpt_scene(p, traj)
display(fig)
readline()
