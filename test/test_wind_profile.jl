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

    @testset "steady_wind" begin
        f = steady_wind(8.5)
        @test f(0.0)   ≈ 8.5
        @test f(60.0)  ≈ 8.5
        @test f(120.0) ≈ 8.5
    end

    @testset "wind_ramp" begin
        f = wind_ramp(6.0, 14.0, 30.0, 90.0)
        # Before ramp
        @test f(0.0)  ≈ 6.0
        @test f(30.0) ≈ 6.0
        # After ramp
        @test f(90.0)  ≈ 14.0
        @test f(120.0) ≈ 14.0
        # Midpoint (linear)
        @test f(60.0) ≈ 10.0  atol=1e-10
        # Monotonically increasing through ramp
        vs = [f(t) for t in 30.0:10.0:90.0]
        @test issorted(vs)
    end

    @testset "gust_event" begin
        f = gust_event(8.0, 15.0, 40.0, 60.0)
        # Outside gust window
        @test f(0.0)  ≈ 8.0
        @test f(39.9) ≈ 8.0
        @test f(60.1) ≈ 8.0
        # Peak at midpoint
        v_peak = f(50.0)
        @test v_peak ≈ 15.0  atol=0.01
        # C¹ continuity at edges (derivative near 0 at t_start and t_end)
        @test f(40.0) ≈ 8.0  atol=1e-6
        @test f(60.0) ≈ 8.0  atol=1e-6
        # All values within [v_base, v_gust]
        @test all(8.0 .<= f(t) .<= 15.0 for t in 40.0:0.1:60.0)
    end

    @testset "turbulent_wind" begin
        v_mean = 10.0
        TI     = 0.15
        f      = turbulent_wind(v_mean, TI, 120.0; rng_seed=42)
        ts     = 0.0:1.0:120.0
        vs     = [f(t) for t in ts]
        # Mean within 15 % of v_mean over 120 s
        @test abs(sum(vs) / length(vs) - v_mean) < 0.15 * v_mean
        # All values positive (clamped to ≥ 0.5)
        @test all(v > 0.0 for v in vs)
        # Standard deviation within [50 %, 200 %] of expected σ
        σ_expected = TI * v_mean
        σ_actual   = sqrt(sum((v - v_mean)^2 for v in vs) / length(vs))
        @test 0.5 * σ_expected < σ_actual < 2.0 * σ_expected
        # Deterministic: same seed → same result
        f2 = turbulent_wind(v_mean, TI, 120.0; rng_seed=42)
        @test f2(50.0) ≈ f(50.0)  atol=1e-10
    end

end
