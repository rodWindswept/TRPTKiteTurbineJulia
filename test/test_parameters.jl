using Test

include("../src/parameters.jl")

@testset "SystemParams — 10 kW preset" begin
    p = params_10kw()

    # Geometry — Framework PDF §5.3
    @test p.rotor_radius ≈ 5.0
    @test p.elevation_angle ≈ π/6 atol=1e-6
    @test p.n_blades == 3

    # TRPT geometry — Design Reasoning Report §5.2
    @test p.tether_length ≈ 30.0        # "For a 30m TRPT"
    @test p.n_lines == 5                 # "5 tethers along the length"
    @test p.n_rings == 14                # rings every 2m → (30/2)−1 = 14
    @test p.tether_diameter ≈ 0.003     # 3mm Dyneema type 01505
    @test p.m_ring ≈ 0.4 atol=0.01      # ~400g per ring (CFRP tubes + connectors)

    # Hub altitude consistency: h_ref = tether_length × sin(elevation_angle)
    @test p.h_ref ≈ p.tether_length * sin(p.elevation_angle) atol=0.1

    # Tether material
    @test p.e_modulus ≈ 100e9

    # TRPT shaft mass sanity: rings + tether ≈ 6.6 kg (DRR §5.2)
    tether_mass = p.n_lines * p.tether_length * 0.006  # 6 g/m Dyneema
    shaft_mass  = p.n_rings * p.m_ring + tether_mass
    @test shaft_mass ≈ 6.6 atol=0.5

    # Mass values — Mass Scaling PDF §"Static Lift Kite Mass Bottleneck"
    total_rotor_mass = p.n_blades * p.m_blade
    @test total_rotor_mass ≈ 11.0 atol=0.1

    # Ground station inertia — Mass Scaling PDF §"Drivetrain Mass and Inertia Matching"
    @test p.i_pto ≈ 0.059 atol=0.001

    # Operating parameters — AeroDyn BEM (Rotor_TRTP_Sizing_Iteration2.xlsx)
    @test p.v_wind_ref ≈ 11.0          # rated hub wind speed
    @test p.cp ≈ 0.22 atol=0.01        # NACA4412 3-blade BEM; 0.222–0.234 across sizes
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
