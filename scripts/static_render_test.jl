# scripts/static_render_test.jl — temporary static geometry verification script
using GLMakie
include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/visualization.jl")

p = params_10kw()
alpha_test = π / 4   # 45° total twist — visually interesting

fig = Figure(size=(900, 700))
ax  = Axis3(fig[1, 1],
            title  = "TRPT 3D Geometry — Static Verification (α_tot = 45°)",
            xlabel = "X (m)", ylabel = "Y (m)", zlabel = "Z (m)")

nodes    = compute_trpt_geometry(p, alpha_test)
n_levels = size(nodes, 1)

# Tether lines (blue)
for j in 1:p.n_lines
    xs = nodes[:, j, 1]; ys = nodes[:, j, 2]; zs = nodes[:, j, 3]
    lines!(ax, xs, ys, zs, color=:blue, linewidth=1.5)
end

# Polygon rings — intermediate levels (black)
for i in 2:(n_levels - 1)
    ring_x = [nodes[i, :, 1]; nodes[i, 1, 1]]
    ring_y = [nodes[i, :, 2]; nodes[i, 1, 2]]
    ring_z = [nodes[i, :, 3]; nodes[i, 1, 3]]
    lines!(ax, ring_x, ring_y, ring_z, color=:black, linewidth=1.2)
end

# Rotor ring (red)
rotor_x = [nodes[end, :, 1]; nodes[end, 1, 1]]
rotor_y = [nodes[end, :, 2]; nodes[end, 1, 2]]
rotor_z = [nodes[end, :, 3]; nodes[end, 1, 3]]
lines!(ax, rotor_x, rotor_y, rotor_z, color=:red, linewidth=3.0)

# Ground anchor
scatter!(ax, [0.0], [0.0], [0.0], color=:green, markersize=20)

save("output/static_render_verification.png", fig)
println("Saved to output/static_render_verification.png")
