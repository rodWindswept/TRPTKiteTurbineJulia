using Test
using LinearAlgebra

include("../src/parameters.jl")
include("../src/geometry.jl")

@testset "TRPT geometry — node positions (inclined shaft)" begin
    p = params_10kw()
    alpha_tot = 0.3

    nodes = compute_trpt_geometry(p, alpha_tot)

    @test size(nodes) == (p.n_rings + 2, p.n_lines, 3)
    # Ground ring centre is at origin — centroid z of the ring ≈ 0
    # (individual nodes are not at z=0 because rings are perpendicular to inclined shaft)
    @test sum(nodes[1, :, 3]) / p.n_lines ≈ 0.0 atol=1e-6

    h_hub = p.tether_length * sin(p.elevation_angle)
    # Top ring centre z ≈ h_hub — centroid z of the ring
    @test sum(nodes[end, :, 3]) / p.n_lines ≈ h_hub atol=1e-6

    x_hub = p.tether_length * cos(p.elevation_angle)
    @test nodes[end, 1, 1] ≈ x_hub atol=p.trpt_hub_radius

    n_levels = p.n_rings + 2
    n_seg    = p.n_rings + 1
    l_seg    = p.tether_length / n_seg
    r_bottom = 2.0 * p.tether_length * p.trpt_rL_ratio / n_seg - p.trpt_hub_radius
    shaft_dir = [cos(p.elevation_angle), 0.0, sin(p.elevation_angle)]

    for i in 1:n_levels
        level_idx  = i - 1
        r_expected = r_bottom + level_idx / (n_levels - 1) * (p.trpt_hub_radius - r_bottom)
        ring_centre = level_idx * l_seg .* shaft_dir
        for j in 1:p.n_lines
            node = nodes[i, j, :]
            dist = sqrt(sum((node .- ring_centre).^2))
            @test dist ≈ r_expected atol=1e-6
        end
    end
end

@testset "Zero twist — lines radially aligned" begin
    p         = params_10kw()
    shaft_dir = [0.0, 0.0, 1.0]   # vertical shaft: ring centres on Z axis
    nodes     = compute_trpt_geometry(p, 0.0, shaft_dir)

    # With zero twist and ring centres on the Z axis, world XY angles match
    phi_ground = atan.(nodes[1, :, 2], nodes[1, :, 1])
    phi_top    = atan.(nodes[end, :, 2], nodes[end, :, 1])
    @test phi_ground ≈ phi_top atol=1e-6
end

@testset "Custom shaft_dir — vertical shaft" begin
    p         = params_10kw()
    shaft_dir = [0.0, 0.0, 1.0]
    nodes     = compute_trpt_geometry(p, 0.0, shaft_dir)

    n_seg  = p.n_rings + 1
    @test all(nodes[1, :, 3] .≈ 0.0)
    @test nodes[end, 1, 3] ≈ p.tether_length atol=1e-6

    n_levels = p.n_rings + 2
    r_bottom = 2.0 * p.tether_length * p.trpt_rL_ratio / n_seg - p.trpt_hub_radius
    for i in 1:n_levels
        level_idx  = i - 1
        r_expected = r_bottom + level_idx / (n_levels - 1) * (p.trpt_hub_radius - r_bottom)
        for j in 1:p.n_lines
            r_actual = sqrt(nodes[i, j, 1]^2 + nodes[i, j, 2]^2)
            @test r_actual ≈ r_expected atol=1e-6
        end
    end
end

@testset "Blade geometry — shape and dimensions" begin
    p         = params_10kw()
    alpha_tot = 0.0
    blades    = compute_blade_geometry(p, alpha_tot)

    @test size(blades) == (p.n_blades, 4, 3)

    BLADE_INNER_FRAC = 0.30
    blade_span    = (p.rotor_radius - p.trpt_hub_radius) / (1.0 - BLADE_INNER_FRAC)
    blade_inner_r = p.trpt_hub_radius - BLADE_INNER_FRAC * blade_span
    blade_outer_r = p.rotor_radius

    n_seg        = p.n_rings + 1
    l_seg        = p.tether_length / n_seg
    shaft_dir    = [cos(p.elevation_angle), 0.0, sin(p.elevation_angle)]
    rotor_centre = (p.n_rings + 1) * l_seg .* shaft_dir

    for b in 1:p.n_blades
        inner_mid = (blades[b, 1, :] .+ blades[b, 2, :]) ./ 2
        outer_mid = (blades[b, 3, :] .+ blades[b, 4, :]) ./ 2
        inner_dist = sqrt(sum((inner_mid .- rotor_centre).^2))
        outer_dist = sqrt(sum((outer_mid .- rotor_centre).^2))
        @test inner_dist ≈ blade_inner_r atol=1e-6
        @test outer_dist ≈ blade_outer_r atol=1e-6
    end
end

@testset "Ground plane — grid structure" begin
    gp = world_ground_plane()
    @test length(gp) > 0
    for (xs, ys, zs) in gp
        @test length(xs) == 2
        @test length(ys) == 2
        @test all(zs .≈ 0.0)
    end
end
