# scripts/record_simulation.jl
# Render the full 120 s TRPT simulation as an MP4 dashboard video.
# Usage: julia --project=. scripts/record_simulation.jl

using GLMakie, DifferentialEquations, Statistics

# Point Makie at the Julia-bundled ffmpeg
using FFMPEG_jll
ENV["FFMPEG"] = FFMPEG_jll.ffmpeg_path

include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/aerodynamics.jl")
include("../src/dynamics.jl")
include("../src/geometry.jl")
include("../src/force_analysis.jl")
include("../src/visualization.jl")

mkpath("output")

# ── 1. Solve ODE ──────────────────────────────────────────────────────────────
println("Solving 120 s ODE...")
p    = params_10kw()
u0   = [0.0, 1.0]
prob = ODEProblem(trpt_ode!, u0, (0.0, 120.0), p)
sol  = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=1/30)
println("  $(length(sol.t)) frames  retcode: $(sol.retcode)")

h_hub = hub_altitude(p.tether_length, p.elevation_angle)
v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h_hub)
n_seg = p.n_rings + 1

traj = (
    t         = sol.t,
    alpha_tot = sol[1, :],
    omega     = sol[2, :],
    power_kw  = [instantaneous_power(p, [sol[1,i], sol[2,i]]) / 1000.0
                 for i in eachindex(sol.t)],
    v_hub     = v_hub,
)

# ── 2. Build scene ────────────────────────────────────────────────────────────
println("Building GLMakie scene...")
fig, time_obs = build_trpt_scene(p, traj)
n_frames = length(traj.t)

# ── 3. Record ─────────────────────────────────────────────────────────────────
out_path = "output/simulation_dashboard.mp4"
println("Recording $n_frames frames → $out_path ...")

record(fig, out_path, 1:n_frames; framerate=30) do frame_idx
    time_obs[] = frame_idx
end

println("Done — $(out_path)")
