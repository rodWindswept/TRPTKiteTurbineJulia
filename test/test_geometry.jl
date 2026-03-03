using Test

include("../src/parameters.jl")
include("../src/visualization.jl")

@testset "TRPT geometry — node positions" begin
    p = params_10kw()
    alpha_tot = 0.3   # small twist

    nodes = compute_trpt_geometry(p, alpha_tot)

    # Shape: (n_rings + 2) levels × n_lines lines × 3 coords
    @test size(nodes) == (p.n_rings + 2, p.n_lines, 3)

    # Ground level (index 1) must be at z = 0
    @test all(nodes[1, :, 3] .≈ 0.0)

    # Top level must be at h_hub altitude
    h_hub = p.tether_length * sin(p.elevation_angle)
    @test nodes[end, 1, 3] ≈ h_hub atol=1e-6

    # Each level must be at its tapered TRPT radius
    n_seg    = p.n_rings + 1
    r_bottom = 2.0 * p.tether_length * p.trpt_rL_ratio / n_seg - p.trpt_hub_radius
    n_levels = p.n_rings + 2
    for i in 1:n_levels
        level_idx = i - 1
        r_expected = r_bottom + level_idx / (n_levels - 1) * (p.trpt_hub_radius - r_bottom)
        for j in 1:p.n_lines
            r_actual = sqrt(nodes[i, j, 1]^2 + nodes[i, j, 2]^2)
            @test r_actual ≈ r_expected atol=1e-6
        end
    end

    # Ground radius ≈ r_bottom, rotor radius ≈ trpt_hub_radius
    r_ground = sqrt(nodes[1, 1, 1]^2 + nodes[1, 1, 2]^2)
    r_top    = sqrt(nodes[end, 1, 1]^2 + nodes[end, 1, 2]^2)
    @test r_ground ≈ r_bottom atol=1e-6
    @test r_top    ≈ p.trpt_hub_radius atol=1e-6
end

@testset "Zero twist — lines radially aligned" begin
    p     = params_10kw()
    nodes = compute_trpt_geometry(p, 0.0)

    # With zero twist all levels have the same angular positions
    phi_ground = atan.(nodes[1, :, 2], nodes[1, :, 1])
    phi_top    = atan.(nodes[end, :, 2], nodes[end, :, 1])
    @test phi_ground ≈ phi_top atol=1e-6
end
