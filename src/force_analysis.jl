# src/force_analysis.jl
# Per-element structural force analysis derived analytically from ODE state.
# All forces in SI units (N, N·m).

const TETHER_SWL = 3500.0   # Dyneema 3mm safe working load (N) — DRR §5.2
const RING_SWL   = 500.0    # CFRP tube ring strut conservative buckling limit (N)

"""
    ForceState

Per-element structural forces at one ODE state.
"""
struct ForceState
    tether_tension   :: Vector{Float64}   # axial tension per segment (N), length n_seg
    ring_compression :: Vector{Float64}   # hoop compression per ring (N), length n_rings
    tau_transmitted  :: Float64           # torque through TRPT (N·m)
    tau_aero         :: Float64           # aerodynamic torque (N·m)
    tau_drag         :: Float64           # tether drag torque (N·m)
end

"""
    element_forces(p, u, v_hub) -> ForceState

Derive per-element structural forces from ODE state `u = [alpha_tot, omega]`.

Tether tension per segment i has four contributions:
  1. Lifter kite pre-tension: vertical component of the lifter line compensates full
     airborne weight. T_lifter = m_airborne × g / (sin(lifter_elevation) × n_lines).
     Same for every segment (top-applied, propagates down all tethers).
  2. Aerodynamic rotor thrust: F_thrust = P_aero / v_eff, distributed to n_lines.
     Same for every segment (top-applied, propagates down all tethers).
  3. Centrifugal: varies per segment — tether mass × ω² × r_i / n_lines.
  4. Gravitational: weight of all mass above this segment along shaft axis / n_lines.

Ring hoop compression from polygon geometry:
  C_i = T_i × r_i / (2 × sin(π / n_lines))
"""
function element_forces(p::SystemParams, u::Vector{Float64}, v_hub::Float64)::ForceState
    n_seg    = p.n_rings + 1
    r_top    = p.trpt_hub_radius
    r_bottom = 2.0 * p.tether_length * p.trpt_rL_ratio / n_seg - r_top
    l_seg    = p.tether_length / n_seg

    spring_coeff = p.e_modulus * π * (p.tether_diameter / 2)^2 * p.n_lines * p.trpt_rL_ratio

    alpha_tot = u[1]
    omega     = u[2]
    g         = 9.81

    tether_line_mass_per_m = p.rho * π * (p.tether_diameter / 2)^2
    m_seg_tether = tether_line_mass_per_m * l_seg * p.n_lines

    tether_tension   = zeros(Float64, n_seg)
    ring_compression = zeros(Float64, p.n_rings)

    # τ_transmitted
    sum_inv_k = sum(1.0 / (spring_coeff * (r_bottom + i / (n_seg - 1) * (r_top - r_bottom)))
                    for i in 0:n_seg-1)
    k_eff = 1.0 / sum_inv_k
    tau_transmitted = k_eff * sin(alpha_tot / n_seg)

    # Aerodynamic torque
    omega_safe = max(abs(omega), 0.1)
    P_aero     = 0.5 * p.rho * v_hub^3 * π * p.rotor_radius^2 * p.cp * cos(p.elevation_angle)^3
    tau_aero   = sign(omega) * P_aero / omega_safe

    # Tether drag torque
    lambda_t   = omega * p.rotor_radius / max(v_hub, 0.1)
    V_a        = v_hub * (lambda_t + sin(p.elevation_angle))
    drag_force = 0.25 * 1.0 * p.tether_diameter * p.tether_length * p.rho * V_a^2
    tau_drag   = drag_force * p.rotor_radius * 0.5

    # ── Top-applied axial loads (same for every segment) ──────────────────────

    # 1. Lifter kite pre-tension: vertical component compensates full airborne weight.
    #    Lifting line at p.lifter_elevation above horizontal; force distributed to n_lines nodes.
    m_airborne    = p.n_blades * p.m_blade + p.n_rings * p.m_ring
    T_lift_total  = m_airborne * g / sin(p.lifter_elevation)
    T_lifter      = T_lift_total / p.n_lines

    # 2. Aerodynamic rotor thrust: applied at rotor end, propagates to all segments.
    #    F_thrust = P_aero / v_effective where v_effective = v_hub * cos(β).
    v_eff    = max(v_hub * cos(p.elevation_angle), 0.1)
    F_thrust = P_aero / v_eff
    T_aero   = F_thrust / p.n_lines

    # ── Per-segment loads (vary with position) ────────────────────────────────

    for i in 0:n_seg-1
        r_i = r_bottom + i / (n_seg - 1) * (r_top - r_bottom)

        # Centrifugal load on this segment's tether mass
        F_centrifugal = m_seg_tether * omega^2 * r_i / p.n_lines

        # Gravitational: weight of all mass above this segment along shaft axis
        n_above       = n_seg - 1 - i
        m_rings_above = n_above * p.m_ring
        m_blades_above = (i == n_seg - 1) ? p.n_blades * p.m_blade : 0.0
        m_tether_above = n_above * m_seg_tether
        m_above        = m_rings_above + m_blades_above + m_tether_above
        F_gravity      = m_above * g * sin(p.elevation_angle) / p.n_lines

        tether_tension[i + 1] = max(0.0, T_lifter + T_aero + F_centrifugal + F_gravity)

        if i < p.n_rings
            ring_compression[i + 1] = tether_tension[i + 1] * r_i /
                                       (2.0 * sin(π / p.n_lines))
        end
    end

    return ForceState(tether_tension, ring_compression, tau_transmitted, tau_aero, tau_drag)
end

"""
    run_force_scan(p, u_frames, v_hub_vec) -> (T_max, C_max)

Scan full trajectory for peak tether tension and ring compression.
Used to calibrate the force colourbar scale before playback.
"""
function run_force_scan(p::SystemParams,
                        u_frames::Vector{Vector{Float64}},
                        v_hub_vec::Vector{Float64})
    T_max = 0.0
    C_max = 0.0
    for (u, v) in zip(u_frames, v_hub_vec)
        fs    = element_forces(p, u, v)
        T_max = max(T_max, maximum(fs.tether_tension))
        C_max = max(C_max, maximum(fs.ring_compression))
    end
    return T_max, C_max
end
