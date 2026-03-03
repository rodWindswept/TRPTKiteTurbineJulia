using Test

include("../src/parameters.jl")

@testset "SystemParams — 10 kW preset" begin
    p = params_10kw()

    # Geometry — Framework PDF §5.3
    @test p.rotor_radius ≈ 5.0
    @test p.tether_length ≈ 150.0
    @test p.n_lines == 6
    @test p.n_rings == 10
    @test p.n_blades == 3
    @test p.elevation_angle ≈ π/6 atol=1e-6

    # Tether material
    @test p.tether_diameter ≈ 0.004
    @test p.e_modulus ≈ 100e9

    # Mass values — Mass Scaling PDF §"Static Lift Kite Mass Bottleneck"
    total_rotor_mass = p.n_blades * p.m_blade
    @test total_rotor_mass ≈ 11.0 atol=0.1

    # Ground station inertia — Mass Scaling PDF §"Drivetrain Mass and Inertia Matching"
    @test p.i_pto ≈ 0.059 atol=0.001

    # Operating parameters
    @test p.cp ≈ 0.15
    @test p.c_pto ≈ 5000.0
end

@testset "SystemParams — 50 kW preset" begin
    p50 = params_50kw()
    p10 = params_10kw()

    # Radius scales geometrically: R ∝ (P_ratio)^(1/3)
    scale_factor = (50.0 / 10.0)^(1/3)
    @test p50.rotor_radius ≈ p10.rotor_radius * scale_factor atol=0.1
    @test p50.tether_length ≈ p10.tether_length * scale_factor atol=1.0

    # Total rotor mass scales at 1.35 exponent
    mass_ratio = (50.0 / 10.0)^1.35
    expected_rotor = 11.0 * mass_ratio
    @test p50.n_blades * p50.m_blade ≈ expected_rotor atol=1.0

    # 50 kW system must be larger than 10 kW
    @test p50.rotor_radius > p10.rotor_radius
    @test p50.tether_length > p10.tether_length
end

@testset "mass_scale function" begin
    p10 = params_10kw()
    p20 = mass_scale(p10, 10.0, 20.0)

    # Mass exponent 1.35
    expected_factor = (20.0 / 10.0)^1.35
    @test p20.m_blade ≈ p10.m_blade * expected_factor atol=0.01

    # Radius scales geometrically
    r_factor = (20.0 / 10.0)^(1/3)
    @test p20.rotor_radius ≈ p10.rotor_radius * r_factor atol=0.01

    # Physics constants should not change
    @test p20.cp == p10.cp
    @test p20.rho == p10.rho
    @test p20.elevation_angle == p10.elevation_angle
end
