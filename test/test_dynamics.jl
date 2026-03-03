using Test
include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/dynamics.jl")

@testset "trpt_ode! — basic call" begin
    p  = params_10kw()
    du = zeros(2)
    u  = [0.1, 1.0]

    trpt_ode!(du, u, p, 0.0)

    # Both derivatives must be finite Float64 values
    @test isfinite(du[1])
    @test isfinite(du[2])
    @test du[1] isa Float64
    @test du[2] isa Float64
end

@testset "instantaneous_power — zero twist" begin
    p = params_10kw()

    # At zero twist α_seg = 0 → sin(0) = 0 → τ_transmitted = 0 → P = 0
    @test instantaneous_power(p, [0.0, 1.0]) ≈ 0.0
end

@testset "instantaneous_power — collapse guard" begin
    p = params_10kw()

    # α_seg = 100.0 / (n_rings + 1) = 100/15 ≈ 6.67 rad ≫ 0.95π ≈ 2.98 rad
    # Collapse guard must force τ_transmitted = 0 → P = 0
    @test instantaneous_power(p, [100.0, 1.0]) ≈ 0.0
end

@testset "instantaneous_power — positive under normal operation" begin
    p = params_10kw()

    # u = [0.5, 2.0] — modest twist, healthy rotation speed
    @test instantaneous_power(p, [0.5, 2.0]) > 0.0
end

@testset "trpt_ode! — rotor accelerates from near-zero ω" begin
    p  = params_10kw()
    du = zeros(2)

    # Very low ω, zero twist → aerodynamic torque dominates → rotor should accelerate
    u = [0.0, 0.01]
    trpt_ode!(du, u, p, 0.0)

    @test du[2] > 0.0
end

@testset "I_total — positive inertia" begin
    p = params_10kw()

    I_total = p.n_blades * p.m_blade * p.rotor_radius^2 +
              p.n_rings  * p.m_ring  * p.rotor_radius^2 +
              p.i_pto

    @test I_total > 0.0
end
