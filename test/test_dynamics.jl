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

    # Replicate the tapered ring inertia calculation from trpt_ode!
    n_seg   = p.n_rings + 1
    r_top_i = p.trpt_hub_radius
    r_bot_i = 2.0 * p.tether_length * p.trpt_rL_ratio / n_seg - r_top_i
    I_rings = sum(p.m_ring * (r_bot_i + i / p.n_rings * (r_top_i - r_bot_i))^2
                  for i in 1:p.n_rings)
    I_total = p.n_blades * p.m_blade * p.rotor_radius^2 + I_rings + p.i_pto

    # Sanity checks
    @test I_rings > 0.0
    @test I_total > 0.0

    # Structural check: I_total is the sum of its three components
    @test I_total ≈ p.n_blades * p.m_blade * p.rotor_radius^2 + I_rings + p.i_pto

    # Blade contribution (3 × 11/3 kg × 5² m²) ≈ 275.0 kg·m²
    @test p.n_blades * p.m_blade * p.rotor_radius^2 ≈ 275.0 atol=1.0

    # Tapered ring inertia must be less than if all rings were at rotor_radius
    @test I_rings < p.n_rings * p.m_ring * p.rotor_radius^2

    # Tapered inertia (rings at 0.96–2.0 m) should be ≈ 13.4 kg·m²
    @test I_rings ≈ 13.4 atol=0.5
end

@testset "TRPT taper — stiffness and inertia" begin
    p = params_10kw()
    n_seg     = p.n_rings + 1
    r_top     = p.trpt_hub_radius
    r_bottom  = 2.0 * p.tether_length * p.trpt_rL_ratio / n_seg - r_top
    coeff     = p.e_modulus * π * (p.tether_diameter / 2)^2 * p.n_lines * p.trpt_rL_ratio

    # Effective stiffness is less than the stiffest (top) segment
    k_top    = coeff * r_top
    k_bottom = coeff * r_bottom
    sum_inv  = sum(1.0 / (coeff * (r_bottom + i / (n_seg - 1) * (r_top - r_bottom)))
                   for i in 0:n_seg-1)
    k_eff    = 1.0 / sum_inv
    @test k_eff < k_top
    @test k_eff < k_bottom   # series combination is always less stiff than weakest element

    # Stress guard reachability: ground-end deflection threshold must be lower than
    # kinematic threshold — i.e. k_bottom/k_eff < n_seg.
    # If this fails the secondary collapse guard is unreachable for this parameter set.
    @test k_bottom / k_eff < n_seg

    # Power is positive at normal operating state
    @test instantaneous_power(p, [0.5, 2.0]) > 0.0

    # ODE call still works and returns finite values
    du = zeros(2)
    trpt_ode!(du, [0.1, 1.0], p, 0.0)
    @test isfinite(du[1])
    @test isfinite(du[2])
end
