using Test
include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/aerodynamics.jl")
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

    # Blade contribution (5 × 11/3 kg × 5² m²) ≈ 458.3 kg·m²
    @test p.n_blades * p.m_blade * p.rotor_radius^2 ≈ 5 * (11.0/3.0) * 25.0 atol=1.0

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

@testset "trpt_ode_limited! — basic call (3-state)" begin
    p  = params_10kw()
    du = zeros(3)
    u  = [0.1, 1.0, deg2rad(23.0)]   # α, ω, β at minimum

    trpt_ode_limited!(du, u, p, 0.0)

    @test isfinite(du[1])
    @test isfinite(du[2])
    @test isfinite(du[3])
end

@testset "trpt_ode_limited! — β stays at floor under low wind" begin
    p  = params_10kw()
    # Low wind (4 m/s) → P_gen < P_rated → limiter tries to decrease β → clamps at β_min
    using DifferentialEquations
    sol = solve(ODEProblem(trpt_ode_limited!,
                           [0.0, 4.1*4.0/p.rotor_radius*0.8, p.β_min],
                           (0.0, 30.0), p),
                Tsit5(); reltol=1e-6, abstol=1e-6, saveat=1.0)
    β_trace = sol[3, :]
    @test all(β_trace .>= p.β_min - 1e-4)   # never goes below floor
end

@testset "trpt_ode_limited! — β rises under high wind" begin
    p  = params_10kw()
    # High wind (18 m/s): P_aero >> P_rated → β must rise above β_min
    using DifferentialEquations
    ω0 = 4.1 * 18.0 / p.rotor_radius
    sol = solve(ODEProblem(trpt_ode_limited!,
                           [0.0, ω0, p.β_min],
                           (0.0, 90.0), p),
                Tsit5(); reltol=1e-6, abstol=1e-6, saveat=1.0)
    β_trace = sol[3, :]
    @test maximum(β_trace) > p.β_min + deg2rad(1.5)   # β rises meaningfully (>1.5° above floor)
    @test all(β_trace .<= p.β_max + 1e-4)              # never exceeds ceiling
end

@testset "trpt_ode_limited! — power clamped near rated at high wind" begin
    p  = params_10kw()
    using DifferentialEquations, Statistics
    ω0 = 4.1 * 18.0 / p.rotor_radius
    sol = solve(ODEProblem(trpt_ode_limited!,
                           [0.0, ω0, p.β_min],
                           (0.0, 120.0), p),
                Tsit5(); reltol=1e-6, abstol=1e-6, saveat=1.0)
    # Average power over final 30 s must be within 15% of rated (limiter settling)
    idx = findall(t -> t >= 90.0, sol.t)
    P_vals = [instantaneous_power(p, [sol[1,i], sol[2,i]], sol[3,i]) for i in idx]
    @test mean(P_vals) <= p.p_rated_w * 1.15
end

@testset "instantaneous_power — 3-arg variant" begin
    p = params_10kw()
    # The 3-arg variant accepts β for type-compatibility with 3-state trajectories.
    # Ground power = τ_transmitted × ω_ground — computed from mechanical state (α, ω),
    # so β does not change the result for the same (α, ω) state.
    P_default = instantaneous_power(p, [0.5, 2.0])
    P_with_β  = instantaneous_power(p, [0.5, 2.0], deg2rad(60.0))
    @test P_default isa Float64
    @test P_with_β  isa Float64
    @test isfinite(P_with_β)
    @test P_with_β == P_default   # β doesn't enter the ground-power formula
end
