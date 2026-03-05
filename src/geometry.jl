# src/geometry.jl
# 3D geometry for the TRPT tensegrity structure.
# Coordinate system: +x downwind, +y crosswind-right, +z vertical up.
# TRPT ground anchor at world origin (0,0,0).
#
# Ring planes are PERPENDICULAR to the shaft axis (orthogonal to shaft_dir).
# Ring centres are spaced along the inclined shaft axis.

using LinearAlgebra

const BLADE_INNER_FRAC = 0.30   # fraction of blade span inboard of the TRPT ring

"""
    compute_trpt_geometry(p, alpha_tot[, shaft_dir])

Compute 3D node positions for the TRPT tensegrity at a given total twist angle.

`shaft_dir` is a unit vector along the TRPT shaft axis (default: derived from
`p.elevation_angle` — `[cos(β), 0, sin(β)]`, pointing downwind).

Ring planes are horizontal (parallel to ground); ring centres are placed along the
shaft axis. This means each node position is:

    nodes[i,j] = ring_centre_i + r_i * [cos(phi_ij), sin(phi_ij), 0]

where `ring_centre_i = (i-1) * l_seg * shaft_dir` and `phi_ij` accumulates twist.

Returns array of size `(n_levels, n_lines, 3)`:
- Axis 1: ring levels, index 1 = ground anchor, index end = rotor end
- Axis 2: tether lines
- Axis 3: [x, y, z] world coordinates (metres)
"""
function compute_trpt_geometry(p::SystemParams, alpha_tot::Float64,
                                shaft_dir::Vector{Float64} = [cos(p.elevation_angle),
                                                              0.0,
                                                              sin(p.elevation_angle)])
    n_levels  = p.n_rings + 2
    n_seg     = p.n_rings + 1
    l_seg     = p.tether_length / n_seg
    alpha_seg = alpha_tot / n_seg

    r_bottom = 2.0 * p.tether_length * p.trpt_rL_ratio / n_seg - p.trpt_hub_radius

    # Perpendicular basis in the plane orthogonal to shaft_dir
    ref   = abs(shaft_dir[3]) < 0.99 ? [0.0, 0.0, 1.0] : [0.0, 1.0, 0.0]
    perp1 = normalize(cross(shaft_dir, ref))
    perp2 = cross(shaft_dir, perp1)

    nodes = zeros(Float64, n_levels, p.n_lines, 3)

    for i in 1:n_levels
        level_idx   = i - 1
        r_i         = r_bottom + level_idx / (n_levels - 1) * (p.trpt_hub_radius - r_bottom)
        ring_centre = level_idx * l_seg .* shaft_dir
        theta_i     = level_idx * alpha_seg

        for j in 1:p.n_lines
            phi        = theta_i + (j - 1) * (2π / p.n_lines)
            node       = ring_centre .+ r_i .* (cos(phi) .* perp1 .+ sin(phi) .* perp2)
            nodes[i, j, :] = node
        end
    end

    return nodes
end

"""
    compute_blade_geometry(p, alpha_tot[, shaft_dir])

Compute 3D blade quad corners for all rotor blades at a given total twist angle.

Returns array of size `(n_blades, 4, 3)` — four corners per blade.

The TRPT ring sits at trpt_hub_radius; ~30% of the blade span is inboard of the ring.
Blades lie in the horizontal rotor plane at hub altitude.
"""
function compute_blade_geometry(p::SystemParams, alpha_tot::Float64,
                                 shaft_dir::Vector{Float64} = [cos(p.elevation_angle),
                                                               0.0,
                                                               sin(p.elevation_angle)])
    blade_span    = (p.rotor_radius - p.trpt_hub_radius) / (1.0 - BLADE_INNER_FRAC)
    blade_inner_r = p.trpt_hub_radius - BLADE_INNER_FRAC * blade_span
    blade_outer_r = p.rotor_radius
    chord         = p.rotor_radius * 0.15

    n_seg        = p.n_rings + 1
    l_seg        = p.tether_length / n_seg
    rotor_centre = (p.n_rings + 1) * l_seg .* shaft_dir

    # Perpendicular basis in the rotor plane (orthogonal to shaft_dir)
    ref   = abs(shaft_dir[3]) < 0.99 ? [0.0, 0.0, 1.0] : [0.0, 1.0, 0.0]
    perp1 = normalize(cross(shaft_dir, ref))
    perp2 = cross(shaft_dir, perp1)

    blades = zeros(Float64, p.n_blades, 4, 3)

    for b in 1:p.n_blades
        phi_b = alpha_tot + (b - 1) * (2π / p.n_blades)
        # blade radial direction — in the plane perpendicular to shaft
        blade_dir = cos(phi_b) .* perp1 .+ sin(phi_b) .* perp2
        # chord direction — perpendicular to blade, within rotor plane
        chord_dir = cross(shaft_dir, blade_dir)

        inner_centre = rotor_centre .+ blade_inner_r .* blade_dir
        outer_centre = rotor_centre .+ blade_outer_r .* blade_dir

        blades[b, 1, :] = inner_centre .+ (chord / 2) .* chord_dir
        blades[b, 2, :] = inner_centre .- (chord / 2) .* chord_dir
        blades[b, 3, :] = outer_centre .- (chord / 2) .* chord_dir
        blades[b, 4, :] = outer_centre .+ (chord / 2) .* chord_dir
    end

    return blades
end

"""
    world_ground_plane(half_extent=40.0, step=5.0)

Return a vector of line segment tuples `(xs, ys, zs)` forming a grid at z = 0.
"""
function world_ground_plane(half_extent::Float64 = 40.0, step::Float64 = 5.0)
    result = Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}[]
    coords = collect(-half_extent:step:half_extent)
    for c in coords
        push!(result, ([-half_extent, half_extent], [c, c],              [0.0, 0.0]))
        push!(result, ([c, c],              [-half_extent, half_extent], [0.0, 0.0]))
    end
    return result
end
