# src/dynamics.jl
# TRPT ODE system for the Kite Turbine Simulator.
#
# Physics references:
#   - Tulloch tether drag model: Tulloch & Yim (2008), AWEC proceedings
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
2. Total rotational inertia of the airborne assembly plus PTO.
3. Spatial discretization of the TRPT shaft into n_rings+1 segments.
4. Torsional stiffness torque transmitted through one segment.
5. Collapse guard: zero transmitted torque when |α_seg| ≥ 0.95π.
6. Aerodynamic power and torque on the rotor disc.
7. Tether drag (Tulloch model), expressed as an equivalent braking torque.
8. Ground PTO angular velocity derived from damping law: τ = c_pto × ω_ground.
9. Equations of motion assembled from the above quantities.
"""
function trpt_ode!(du, u, p::SystemParams, t)
    # Step 1 — Wind speed at hub altitude
    h     = hub_altitude(p.tether_length, p.elevation_angle)
    v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h)

    # Step 2 — Total rotational inertia (blades + rings + PTO)
    I_total = p.n_blades * p.m_blade * p.rotor_radius^2 +
              p.n_rings  * p.m_ring  * p.rotor_radius^2 +
              p.i_pto

    # Step 3 — TRPT spatial discretization
    n_segments = p.n_rings + 1
    α_seg      = u[1] / n_segments        # twist per segment (rad)
    L_seg      = p.tether_length / n_segments

    # Step 4 — Segmented torsional stiffness torque
    k_seg         = (p.e_modulus * π * (p.tether_diameter / 2)^2 *
                     p.n_lines * p.rotor_radius^2) / L_seg
    τ_transmitted = k_seg * sin(α_seg)

    # Step 5 — Collapse guard (tensegrity structural collapse at |α_seg| ≥ 0.95π)
    if abs(α_seg) >= 0.95 * π
        τ_transmitted = 0.0
    end

    # Step 6 — Aerodynamic torque
    ω      = u[2]
    ω_safe = max(abs(ω), 0.1)   # guard against division by zero
    P_aero = 0.5 * p.rho * v_hub^3 * π * p.rotor_radius^2 *
             p.cp * cos(p.elevation_angle)^3
    τ_aero = P_aero / ω_safe

    # Step 7 — Tether drag (Tulloch model)
    λ_t       = ω * p.rotor_radius / max(v_hub, 0.1)   # tether speed ratio
    V_a       = v_hub * (λ_t + sin(p.elevation_angle))  # apparent velocity (m/s)
    drag_force = 0.25 * 1.0 * p.tether_diameter * p.tether_length * p.rho * V_a^2
    τ_drag    = drag_force * p.rotor_radius * 0.5

    # Step 8 — Ground PTO angular velocity (damping law: τ = c_pto × ω_ground)
    ω_ground = τ_transmitted / p.c_pto

    # Step 9 — Equations of motion
    du[1] = ω - ω_ground                                          # dα_tot/dt
    du[2] = (τ_aero - τ_drag - τ_transmitted) / I_total          # dω/dt

    return nothing
end

"""
    instantaneous_power(p::SystemParams, u::Vector{Float64}) -> Float64

Return the electrical power extracted by the PTO (W) at ODE state `u`.

Computed as P = τ_transmitted × ω_ground, where:
- τ_transmitted is the torsional torque passing through one TRPT segment,
- ω_ground = τ_transmitted / c_pto is the PTO shaft angular velocity.

Applies the collapse guard: returns 0 W when |α_seg| ≥ 0.95π.
"""
function instantaneous_power(p::SystemParams, u::Vector{Float64})::Float64
    n_segments    = p.n_rings + 1
    α_seg         = u[1] / n_segments
    L_seg         = p.tether_length / n_segments
    k_seg         = (p.e_modulus * π * (p.tether_diameter / 2)^2 *
                     p.n_lines * p.rotor_radius^2) / L_seg
    τ_transmitted = k_seg * sin(α_seg)

    if abs(α_seg) >= 0.95 * π
        τ_transmitted = 0.0
    end

    ω_ground = τ_transmitted / p.c_pto
    return τ_transmitted * ω_ground
end
