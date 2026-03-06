using Test

include("../src/parameters.jl")

@testset "SystemParams — 10 kW preset" begin
    p = params_10kw()

    # Geometry — Framework PDF §5.3
    @test p.rotor_radius ≈ 5.0
    @test p.elevation_angle ≈ π/6 atol=1e-6
    @test p.n_blades == 5   # one blade per tether line / polygon vertex

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

    # Mass values — 5 blades at 11/3 kg each (DRR blade design mass, 5 blades from n_lines=5)
    total_rotor_mass = p.n_blades * p.m_blade
    @test total_rotor_mass ≈ 5 * (11.0 / 3.0) atol=0.1

    # Ground station inertia — Mass Scaling PDF §"Drivetrain Mass and Inertia Matching"
    @test p.i_pto ≈ 0.059 atol=0.001

    # Operating parameters — AeroDyn BEM (Rotor_TRTP_Sizing_Iteration2.xlsx)
    @test p.v_wind_ref ≈ 11.0          # rated hub wind speed
    @test p.cp ≈ 0.22 atol=0.01        # NACA4412 3-blade BEM; 0.222–0.234 across sizes
    @test p.c_pto ≈ 100.0   # tuned for optimal TSR λ≈4.1 at rated wind (was 5000 pre-BEM)
end

@testset "SystemParams — 50 kW preset" begin
    p50 = params_50kw()
    p10 = params_10kw()

    # Radius scales aerodynamically: R ∝ (P_ratio)^(1/2) since P_aero ∝ R²
    scale_factor = (50.0 / 10.0)^(1/2)
    @test p50.rotor_radius ≈ p10.rotor_radius * scale_factor atol=0.1
    @test p50.tether_length ≈ p10.tether_length * scale_factor atol=1.0

    # Total rotor mass scales at 1.35 exponent (base = 5 blades × 11/3 kg each)
    mass_ratio = (50.0 / 10.0)^1.35
    expected_rotor = (5 * 11.0 / 3.0) * mass_ratio
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

    # Radius scales aerodynamically: R ∝ P^(1/2)
    r_factor = (20.0 / 10.0)^(1/2)
    @test p20.rotor_radius ≈ p10.rotor_radius * r_factor atol=0.01

    # Physics constants should not change
    @test p20.cp == p10.cp
    @test p20.rho == p10.rho
    @test p20.elevation_angle == p10.elevation_angle
end

@testset "TRPT taper geometry — 10 kW preset" begin
    p = params_10kw()

    # rL ratio from DRR Grasshopper
    @test p.trpt_rL_ratio ≈ 0.74 atol=0.01

    # Hub radius at rotor end
    @test p.trpt_hub_radius ≈ 2.0 atol=0.01

    # Derived r_bottom must be positive
    n_seg = p.n_rings + 1
    r_bottom = 2.0 * p.tether_length * p.trpt_rL_ratio / n_seg - p.trpt_hub_radius
    @test r_bottom > 0.0

    # Average segment length = average_r / rL_ratio should ≈ tether_length / n_seg
    avg_r = (p.trpt_hub_radius + r_bottom) / 2.0
    @test avg_r / p.trpt_rL_ratio ≈ p.tether_length / n_seg atol=0.01

    # Mass scaling preserves rL_ratio, scales hub radius as R ∝ P^(1/2)
    p50 = params_50kw()
    @test p50.trpt_rL_ratio ≈ p.trpt_rL_ratio
    scale = (50.0 / 10.0)^(1.0/2.0)
    @test p50.trpt_hub_radius ≈ p.trpt_hub_radius * scale atol=0.01
end
