# scripts/record_50kw.jl
# Record the 50 kW simulation as an MP4 dashboard video.
#
# Usage:
#   julia --project=. scripts/record_50kw.jl

using GLMakie, DifferentialEquations, Statistics, FFMPEG_jll

ENV["FFMPEG"] = FFMPEG_jll.ffmpeg_path

include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/dynamics.jl")
include("../src/geometry.jl")
include("../src/force_analysis.jl")
include("../src/visualization.jl")

mkpath("output")

println("Solving 120 s ODE (50 kW)...")
p   = params_50kw()
sol = solve(ODEProblem(trpt_ode!, [0.0, 1.0], (0.0, 120.0), p),
            Tsit5(); reltol=1e-6, abstol=1e-6, saveat=1/30)
println("  $(length(sol.t)) frames  retcode: $(sol.retcode)")

h_hub = hub_altitude(p.tether_length, p.elevation_angle)
v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h_hub)

traj = (
    t         = sol.t,
    alpha_tot = sol[1, :],
    omega     = sol[2, :],
    power_kw  = [instantaneous_power(p, [sol[1,i], sol[2,i]]) / 1000.0
                 for i in eachindex(sol.t)],
    v_hub     = v_hub,
)

println("Mean $(round(mean(traj.power_kw), digits=2)) kW  |  Peak $(round(maximum(traj.power_kw), digits=2)) kW")

println("Building GLMakie scene...")
fig, time_obs = build_trpt_scene(p, traj)

out_path = "output/simulation_50kw.mp4"
println("Recording $(length(traj.t)) frames → $out_path ...")

record(fig, out_path, 1:length(traj.t); framerate=30) do frame_idx
    time_obs[] = frame_idx
end

println("Done — $out_path")
