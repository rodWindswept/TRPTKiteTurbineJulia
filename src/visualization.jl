# src/visualization.jl
# GLMakie 3D visualization of the TRPT tensegrity structure.
# Strategy: pre-compute ODE trajectory, then animate at ≥30 FPS.

using GLMakie

# ── Geometry ──────────────────────────────────────────────────────────────────

"""
    compute_trpt_geometry(p, alpha_tot)

Compute 3D node positions for the TRPT tensegrity at a given total twist angle.

Returns a 3D array of size `(n_levels, n_lines, 3)` where:
- Axis 1: ring levels from ground (index 1) to rotor (index n_rings+2)
- Axis 2: individual tether lines
- Axis 3: [x, y, z] coordinates in metres

Radii use the TRPT taper: r_bottom at ground end (index 1), trpt_hub_radius at
rotor end (index n_rings+2), linearly interpolated for intermediate rings.

The z-axis is vertical (altitude). Each level is rotated by α_seg relative
to the one below, building up the cumulative twist from ground to rotor.
"""
function compute_trpt_geometry(p::SystemParams, alpha_tot::Float64)
    n_levels  = p.n_rings + 2
    n_seg     = p.n_rings + 1
    l_seg     = p.tether_length / n_seg
    alpha_seg = alpha_tot / n_seg

    r_bottom = 2.0 * p.tether_length * p.trpt_rL_ratio / n_seg - p.trpt_hub_radius

    nodes = zeros(Float64, n_levels, p.n_lines, 3)

    for i in 1:n_levels
        level_idx = i - 1          # 0 at ground, n_rings+1 at rotor
        r_i       = r_bottom + level_idx / (n_levels - 1) * (p.trpt_hub_radius - r_bottom)
        z_pos     = level_idx * l_seg * sin(p.elevation_angle)
        theta_i   = level_idx * alpha_seg   # cumulative twist at this level

        for j in 1:p.n_lines
            phi = theta_i + (j - 1) * (2π / p.n_lines)   # angular position of line j
            nodes[i, j, 1] = r_i * cos(phi)
            nodes[i, j, 2] = r_i * sin(phi)
            nodes[i, j, 3] = z_pos
        end
    end

    return nodes
end
