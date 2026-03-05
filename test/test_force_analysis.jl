using Test

include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/dynamics.jl")
include("../src/force_analysis.jl")

@testset "ForceState — structure" begin
    p     = params_10kw()
    n_seg = p.n_rings + 1
    h_hub = hub_altitude(p.tether_length, p.elevation_angle)
    v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h_hub)
    u     = [0.5, 1.8]

    fs = element_forces(p, u, v_hub)

    @test fs isa ForceState
    @test length(fs.tether_tension)   == n_seg
    @test length(fs.ring_compression) == p.n_rings
    @test fs.tau_transmitted          >= 0.0
    @test fs.tau_aero                 >  0.0
    @test fs.tau_drag                 >= 0.0
end

@testset "ForceState — physical bounds" begin
    p     = params_10kw()
    h_hub = hub_altitude(p.tether_length, p.elevation_angle)
    v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h_hub)
    u     = [0.5, 1.8]

    fs = element_forces(p, u, v_hub)

    @test all(fs.tether_tension   .>= 0.0)
    @test all(fs.ring_compression .>= 0.0)
    @test all(fs.tether_tension   .< 3500.0)   # well under Dyneema 3mm SWL
    # At rated: T ≈ 227 N/line → C ≈ 227 × r_top / (2 × sin(π/5)) ≈ 384 N < 500 N ✓
    @test all(fs.ring_compression .< 500.0)
end

@testset "ForceState — lifter pre-tension at zero wind" begin
    p  = params_10kw()
    u0 = [0.0, 0.0]
    fs = element_forces(p, u0, 0.1)   # near-zero wind, no rotation

    # Lifter kite provides pre-tension even at rest
    m_airborne   = p.n_blades * p.m_blade + p.n_rings * p.m_ring
    T_lift_total = m_airborne * 9.81 / sin(p.lifter_elevation)
    T_expected   = T_lift_total / p.n_lines   # per line, roughly

    @test all(fs.tether_tension .>= T_expected * 0.9)   # within 10% (gravity component varies)
end

@testset "ForceState — zero state" begin
    p     = params_10kw()
    h_hub = hub_altitude(p.tether_length, p.elevation_angle)
    v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h_hub)
    u0    = [0.0, 0.0]

    fs0 = element_forces(p, u0, v_hub)

    @test fs0.tau_transmitted ≈ 0.0 atol=1e-6
    @test all(fs0.tether_tension .>= 0.0)
end

@testset "run_force_scan — max over trajectory" begin
    p     = params_10kw()
    h_hub = hub_altitude(p.tether_length, p.elevation_angle)
    v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h_hub)

    u_frames  = [[0.1, 1.0], [0.5, 1.8], [0.3, 1.5]]
    v_hub_vec = fill(v_hub, 3)

    T_max, C_max = run_force_scan(p, u_frames, v_hub_vec)

    @test T_max > 0.0
    @test C_max > 0.0
    for (u, v) in zip(u_frames, v_hub_vec)
        fs = element_forces(p, u, v)
        @test T_max >= maximum(fs.tether_tension)
        @test C_max >= maximum(fs.ring_compression)
    end
end
