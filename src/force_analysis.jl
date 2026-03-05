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

Tether tension per segment i:
  T_i = torsional + centrifugal + gravitational components

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

    alpha_per_seg = alpha_tot / n_seg

    tether_tension   = zeros(Float64, n_seg)
    ring_compression = zeros(Float64, p.n_rings)

    sum_inv_k = sum(1.0 / (spring_coeff * (r_bottom + i / (n_seg - 1) * (r_top - r_bottom)))
                    for i in 0:n_seg-1)
    k_eff = 1.0 / sum_inv_k
    tau_transmitted = k_eff * sin(alpha_tot / n_seg)

    omega_safe = max(abs(omega), 0.1)
    P_aero     = 0.5 * p.rho * v_hub^3 * π * p.rotor_radius^2 * p.cp * cos(p.elevation_angle)^3
    tau_aero   = sign(omega) * P_aero / omega_safe

    lambda_t   = omega * p.rotor_radius / max(v_hub, 0.1)
    V_a        = v_hub * (lambda_t + sin(p.elevation_angle))
    drag_force = 0.25 * 1.0 * p.tether_diameter * p.tether_length * p.rho * V_a^2
    tau_drag   = drag_force * p.rotor_radius * 0.5

    for i in 0:n_seg-1
        r_i = r_bottom + i / (n_seg - 1) * (r_top - r_bottom)

        # Per-line tether tension: axial structural loads only (N per tether line).
        # Torsional shear is carried by the twisted helix geometry and tracked separately
        # via tau_transmitted; it does not contribute to axial (pull-out) tether tension
        # in the same way as gravity and centrifugal loads.

        # Centrifugal: total centrifugal force of this segment's tether mass divided by n_lines
        F_centrifugal = m_seg_tether * omega^2 * r_i / p.n_lines

        # Gravitational: weight of all mass above this segment projected along shaft axis,
        # divided equally among n_lines
        n_above        = n_seg - 1 - i
        m_rings_above  = n_above * p.m_ring
        m_blades_above = (i == n_seg - 1) ? p.n_blades * p.m_blade : 0.0
        m_tether_above = n_above * m_seg_tether
        m_above        = m_rings_above + m_blades_above + m_tether_above
        F_gravity      = m_above * g * sin(p.elevation_angle) / p.n_lines

        # Per-line tether tension from axial loads (N)
        tether_tension[i + 1] = max(0.0, F_centrifugal + F_gravity)

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
