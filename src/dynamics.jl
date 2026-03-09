# src/dynamics.jl
# TRPT ODE system for the Kite Turbine Simulator.
#
# Physics references:
#   - Tulloch tether drag model: Tulloch & Yue (2008), AWEC proceedings
#   - TRPT stiffness formulation: "Rotary AWES Julia Simulation Framework.pdf" §3–§5
#   - Segmented TRPT discretization: Framework PDF §3, Figure 3

"""
    trpt_ode!(du, u, p::SystemParams, t)

In-place ODE function for the two-state TRPT system.

State vector `u`:
- `u[1]` = α_tot : total TRPT twist angle (rad)
- `u[2]` = ω     : airborne rotor angular velocity (rad/s)

Derivative vector `du`:
- `du[1]` = dα_tot/dt = ω − ω_ground   (rad/s)
- `du[2]` = dω/dt = (τ_aero − τ_drag − τ_transmitted) / I_total   (rad/s²)

Physics pipeline:
1. Wind speed at hub altitude (Hellmann profile).
2. Total rotational inertia: blades (R²) + rings at tapered positions + PTO.
3. TRPT taper geometry: r_top = trpt_hub_radius; r_bottom from tether_length constraint.
4. Effective series torsional stiffness k_eff from tapered segments (k_i ∝ r_i × rL_ratio).
5. Collapse guard: zero torque when average segment twist ≥ 0.95π (kinematic) or
   ground-end segment stress twist ≥ 0.95π (stress).
6. Aerodynamic power and torque on the rotor disc.
7. Tether drag (Tulloch model), expressed as an equivalent braking torque.
8. Ground PTO angular velocity from MPPT law: τ = k_mppt × ω_ground².
9. Equations of motion assembled from the above quantities.
"""
function trpt_ode!(du, u, p::SystemParams, t)
    # Step 1 — Wind speed at hub altitude
    h     = hub_altitude(p.tether_length, p.elevation_angle)
    v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h)

    # Step 2 — Total rotational inertia (blades + rings at taper positions + PTO)
    # Rings at their actual taper positions (ground=0, rotor=n_rings-1)
    # trpt_hub_radius is rotor end; r_bottom computed from taper constraint
    n_segments_inertia = p.n_rings + 1
    r_top_inertia  = p.trpt_hub_radius
    r_bot_inertia  = 2.0 * p.tether_length * p.trpt_rL_ratio / n_segments_inertia - r_top_inertia
    # n_rings intermediate rings (exclude ground and rotor endpoints, index 0 is bottom)
    I_rings = sum(p.m_ring * (r_bot_inertia + i / p.n_rings * (r_top_inertia - r_bot_inertia))^2
                  for i in 1:p.n_rings)
    I_total = p.n_blades * p.m_blade * p.rotor_radius^2 + I_rings + p.i_pto

    # Step 3 — TRPT taper geometry
    n_segments = p.n_rings + 1
    r_top   = p.trpt_hub_radius
    r_bottom = 2.0 * p.tether_length * p.trpt_rL_ratio / n_segments - r_top

    # Step 4 — Effective series torsional stiffness (tapered segments in series)
    # k_i = E × π × (d/2)² × n_lines × r_i × rL_ratio  (since L_i = r_i/rL)
    spring_coeff = p.e_modulus * π * (p.tether_diameter / 2)^2 * p.n_lines * p.trpt_rL_ratio
    # 1/k_eff = Σ(1/k_i) for linear taper: r_i = r_bottom + i/(n_segments-1)×(r_top−r_bottom)
    sum_inv_k = sum(1.0 / (spring_coeff * (r_bottom + i / (n_segments - 1) * (r_top - r_bottom)))
                    for i in 0:n_segments-1)
    k_eff = 1.0 / sum_inv_k

    # Representative twist per segment for sin model
    α_avg = u[1] / n_segments
    τ_transmitted = k_eff * sin(α_avg)

    # Step 5 — Collapse guard: zero torque when structure is over-twisted
    # Primary kinematic guard: average segment twist ≥ 0.95π means segments have wound
    # past the tensegrity stability limit (shaft has collapsed kinematically).
    # Secondary stress guard: ground-end segment (smallest r, softest spring) experiences
    # the most deflection; ground-end deflection = α_total × (k_eff/k_bottom). If that
    # exceeds 0.95π, collapse.
    k_bottom = spring_coeff * r_bottom
    α_bottom = abs(u[1]) * k_eff / k_bottom   # ground-end segment deflection (linear model)
    if abs(α_avg) >= 0.95 * π || α_bottom >= 0.95 * π
        τ_transmitted = 0.0
    end

    # Step 6 — Aerodynamic torque (TSR-dependent Cp from BEM table)
    ω      = u[2]
    ω_safe = max(abs(ω), 0.1)   # guard against ω→0; 0.1 rad/s ≈ 1 RPM
    λ_t    = ω * p.rotor_radius / max(v_hub, 0.1)   # tip speed ratio
    P_aero = 0.5 * p.rho * v_hub^3 * π * p.rotor_radius^2 *
             cp_at_tsr(λ_t) * cos(p.elevation_angle)^3
    # sign(ω) preserves direction: positive ω → accelerating torque,
    # negative ω (solver transient) → braking torque
    τ_aero = sign(ω) * P_aero / ω_safe

    # Step 7 — Tether drag (Tulloch model)
    # V_a = v_hub × (λ_t + sin(β)): sin(β) is the wind component perpendicular to the
    # inclined tether axis (drag-relevant direction); per Framework PDF §3 / design doc §3.3
    V_a       = v_hub * (λ_t + sin(p.elevation_angle))  # apparent velocity (m/s)
    drag_force = 0.25 * 1.0 * p.tether_diameter * p.tether_length * p.rho * V_a^2
    τ_drag    = drag_force * p.rotor_radius * 0.5

    # Step 8 — Ground PTO angular velocity (MPPT law: τ = k_mppt × ω_ground²)
    ω_ground = sqrt(max(τ_transmitted, 0.0) / p.k_mppt)

    # Step 9 — Equations of motion
    du[1] = ω - ω_ground                                          # dα_tot/dt
    du[2] = (τ_aero - τ_drag - τ_transmitted) / I_total          # dω/dt

    return nothing
end

"""
    trpt_ode_wind!(du, u, pw, t)

Variable-wind variant of `trpt_ode!`. Accepts `pw = (p::SystemParams, wind_fn)`
where `wind_fn(t::Float64) -> Float64` returns the reference wind speed at time
`t` (m/s at altitude `p.h_ref`). Physics are identical to `trpt_ode!`.

Use with `ODEProblem(trpt_ode_wind!, u0, tspan, (p, wind_fn))`.
"""
function trpt_ode_wind!(du, u, pw, t)
    p, wind_fn = pw

    # Step 1 — Wind speed at hub altitude (time-varying)
    h       = hub_altitude(p.tether_length, p.elevation_angle)
    v_ref_t = wind_fn(t)
    v_hub   = wind_at_altitude(v_ref_t, p.h_ref, h)

    # Steps 2–9 identical to trpt_ode! ─────────────────────────────────────────
    n_segments_inertia = p.n_rings + 1
    r_top_inertia  = p.trpt_hub_radius
    r_bot_inertia  = 2.0 * p.tether_length * p.trpt_rL_ratio / n_segments_inertia - r_top_inertia
    I_rings = sum(p.m_ring * (r_bot_inertia + i / p.n_rings * (r_top_inertia - r_bot_inertia))^2
                  for i in 1:p.n_rings)
    I_total = p.n_blades * p.m_blade * p.rotor_radius^2 + I_rings + p.i_pto

    n_segments = p.n_rings + 1
    r_top      = p.trpt_hub_radius
    r_bottom   = 2.0 * p.tether_length * p.trpt_rL_ratio / n_segments - r_top

    spring_coeff = p.e_modulus * π * (p.tether_diameter / 2)^2 * p.n_lines * p.trpt_rL_ratio
    sum_inv_k    = sum(1.0 / (spring_coeff * (r_bottom + i / (n_segments - 1) * (r_top - r_bottom)))
                       for i in 0:n_segments-1)
    k_eff = 1.0 / sum_inv_k

    α_avg         = u[1] / n_segments
    τ_transmitted = k_eff * sin(α_avg)

    k_bottom = spring_coeff * r_bottom
    α_bottom = abs(u[1]) * k_eff / k_bottom
    if abs(α_avg) >= 0.95 * π || α_bottom >= 0.95 * π
        τ_transmitted = 0.0
    end

    ω      = u[2]
    ω_safe = max(abs(ω), 0.1)
    λ_t    = ω * p.rotor_radius / max(v_hub, 0.1)
    P_aero = 0.5 * p.rho * v_hub^3 * π * p.rotor_radius^2 *
             cp_at_tsr(λ_t) * cos(p.elevation_angle)^3
    τ_aero = sign(ω) * P_aero / ω_safe

    V_a        = v_hub * (λ_t + sin(p.elevation_angle))
    drag_force = 0.25 * 1.0 * p.tether_diameter * p.tether_length * p.rho * V_a^2
    τ_drag     = drag_force * p.rotor_radius * 0.5

    ω_ground = sqrt(max(τ_transmitted, 0.0) / p.k_mppt)

    du[1] = ω - ω_ground
    du[2] = (τ_aero - τ_drag - τ_transmitted) / I_total

    return nothing
end

"""
    instantaneous_power(p::SystemParams, u::Vector{Float64}) -> Float64

Return the electrical power extracted by the PTO (W) at ODE state `u`.

Computed as P = τ_transmitted × ω_ground, where:
- τ_transmitted is the torsional torque from the effective series stiffness of the
  tapered TRPT (k_eff from segments in series, each with k_i = E π (d/2)² n_lines r_i rL),
- ω_ground = √(τ_transmitted / k_mppt) is the PTO shaft angular velocity.

Applies the collapse guard: returns 0 W when average segment twist ≥ 0.95π (kinematic
collapse) or ground-end segment stress twist ≥ 0.95π (stress-driven collapse).
"""
function instantaneous_power(p::SystemParams, u::Vector{Float64})::Float64
    n_segments   = p.n_rings + 1
    r_top        = p.trpt_hub_radius
    r_bottom     = 2.0 * p.tether_length * p.trpt_rL_ratio / n_segments - r_top
    spring_coeff = p.e_modulus * π * (p.tether_diameter / 2)^2 * p.n_lines * p.trpt_rL_ratio
    sum_inv_k    = sum(1.0 / (spring_coeff * (r_bottom + i / (n_segments - 1) * (r_top - r_bottom)))
                       for i in 0:n_segments-1)
    k_eff        = 1.0 / sum_inv_k
    α_avg        = u[1] / n_segments
    τ_transmitted = k_eff * sin(α_avg)

    k_bottom = spring_coeff * r_bottom
    α_bottom = abs(u[1]) * k_eff / k_bottom   # ground-end segment deflection (linear model)
    if abs(α_avg) >= 0.95 * π || α_bottom >= 0.95 * π
        τ_transmitted = 0.0
    end

    ω_ground = sqrt(max(τ_transmitted, 0.0) / p.k_mppt)
    return τ_transmitted * ω_ground
end
