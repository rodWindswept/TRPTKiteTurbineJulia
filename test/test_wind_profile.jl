using Test

include("../src/parameters.jl")
include("../src/wind_profile.jl")

@testset "wind_profile" begin

    @testset "wind_at_altitude" begin

        # 1. Reference height returns reference speed
        @test wind_at_altitude(10.0, 100.0, 100.0) ≈ 10.0

        # 2. Higher altitude gives more wind
        @test wind_at_altitude(10.0, 100.0, 200.0) > 10.0

        # 3. Lower altitude gives less wind
        @test wind_at_altitude(10.0, 100.0, 50.0) < 10.0

        # 4. Hellmann exponent (1/7) applies correctly
        @test wind_at_altitude(10.0, 100.0, 200.0) ≈ 10.0 * (200.0 / 100.0)^(1.0 / 7.0)  atol=1e-9

        # 5. Non-positive altitudes return 0.0
        @test wind_at_altitude(10.0, 100.0, 0.0)  ≈ 0.0
        @test wind_at_altitude(10.0, 100.0, -5.0) ≈ 0.0

        # 6. Custom Hellmann exponent is respected
        @test wind_at_altitude(10.0, 100.0, 200.0; hellmann_exponent=0.2) ≈ 10.0 * (200.0 / 100.0)^0.2  atol=1e-9

    end

    @testset "hub_altitude" begin

        # 7. Basic geometry: 30 m tether at 30° elevation → 15 m hub altitude
        @test hub_altitude(30.0, π / 6) ≈ 15.0  atol=0.01

        # 8. Integration with params_10kw: hub altitude is 15 m, wind at that altitude is positive
        p = params_10kw()
        h_hub = hub_altitude(p.tether_length, p.elevation_angle)
        @test h_hub ≈ 15.0  atol=0.01
        @test wind_at_altitude(p.v_wind_ref, p.h_ref, h_hub) > 0.0

    end

end
