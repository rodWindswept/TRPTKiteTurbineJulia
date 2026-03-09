# src/parameters.jl
# Physical constants and preset configurations for the TRPT Kite Turbine Simulator.
# All numerical values derived from:
#   - "Rotary AWES Julia Simulation Framework.pdf"         (Framework PDF)
#   - "Kite Turbine Mass Scaling Analysis.pdf"              (Mass Scaling PDF)
#   - "Design Reasoning Report - General Release.pdf" §5.2  (DRR)
#   - data_collection_v3.xlsx                               (rotor geometry survey)
#   - Rotor_TRTP_Sizing_Iteration2.xlsx                     (AeroDyn BEM simulations)

"""
    SystemParams

All physical parameters for one TRPT configuration.
All fields use SI units: m, kg, Pa, rad, rad/s, N·m.
"""
struct SystemParams
    # Atmospheric
    rho::Float64              # Air density (kg/m³), standard 1.225
    v_wind_ref::Float64       # Reference wind speed at h_ref (m/s)
    h_ref::Float64            # Reference altitude for Hellmann wind profile (m)

    # TRPT geometry — Framework PDF §3
    elevation_angle::Float64  # TRPT shaft elevation angle β (rad)
    lifter_elevation::Float64 # Lifter kite line elevation angle (rad); vertical component compensates airborne weight
    rotor_radius::Float64     # Airborne ring radius R (m)
    tether_length::Float64    # Unstretched total TRPT length L₀ (m)

    # TRPT tapered geometry
    trpt_hub_radius::Float64  # Ring radius at the rotor end of the TRPT (m); r_top of taper
    trpt_rL_ratio::Float64    # Geometry ratio r/L per segment (dimensionless).
                              # Each segment: L_i = r_i / trpt_rL_ratio
                              # DRR Grasshopper FEA: rL = 0.740741; Iteration2 sizing: 0.5

    # Tether material — Framework PDF §3
    n_lines::Int64            # Number of TRPT tether lines
    tether_diameter::Float64  # Tether line diameter d_t (m)
    e_modulus::Float64        # Young's modulus (Pa); Dyneema ≈ 100 GPa

    # Polygon spacer rings — Framework PDF §3
    n_rings::Int64            # Number of polygon spacer rings
    m_ring::Float64           # Mass per ring (kg)

    # Lifting blades/kites — Framework PDF §3
    n_blades::Int64           # Number of lifting blades
    m_blade::Float64          # Mass per blade (kg)

    # Aerodynamics
    # rotor_radius is the aerodynamic outer radius used in the power formula
    # (P = ½ρv³πR²Cp cos³β), blade inertia (I = n_blades × m_blade × R²),
    # and tether speed ratio (λ_t = ω × R / v_hub).
    # TRPT stiffness and ring inertia use trpt_hub_radius / trpt_rL_ratio instead.
    cp::Float64               # Rotor power coefficient; AeroDyn BEM ≈ 0.22 (NACA4412, 3-blade)

    # Ground station — Mass Scaling PDF §"Drivetrain Mass and Inertia Matching"
    i_pto::Float64            # Total PTO rotational inertia (kg·m²)
    k_mppt::Float64           # MPPT gain constant k (N·m·s²/rad²): τ_gen = k × ω_ground²
                              # Sets the quadratic load curve for maximum power point tracking.
                              # Derived from rated torque and speed: k = τ_rated / ω_rated²
                              # Eliminates the bistability of a fixed linear damper below rated wind.

    # Power limiter — elevation-angle proportional controller
    p_rated_w  ::Float64   # Rated electrical power (W). Default: 10_000.0
    β_min      ::Float64   # Minimum safe elevation angle (rad). Default: deg2rad(23.0)
                           #   Floor: blade tips must clear the ground at lowest operating angle
    β_max      ::Float64   # Maximum safe elevation angle (rad). Default: deg2rad(67.0)
                           #   Ceiling: lifter kite pull limit (varies IRL with wind speed)
    β_rate_max ::Float64   # Maximum elevation change rate (rad/s). Default: deg2rad(1.0)
                           #   1 °/s — representative of lifter kite / mooring servo speed
    kp_elev    ::Float64   # Proportional gain (rad/W/s). Default: 5e-5
                           #   Gives ≈ 0.5 °/s elevation response per 1 kW overpower
end

"""
    params_10kw()

10 kW prototype configuration — empirically validated from PDFs, DRR, and AeroDyn BEM.

Rated operating point (DRR; AeroDyn sizing data):
  - Rated wind speed : 11 m/s (all AeroDyn simulations in Rotor_TRTP_Sizing_Iteration2.xlsx)
  - Elevation angle  : 30° (as physically built; newer designs optimise at 20°)

Mass budget (Mass Scaling PDF §"Static Lift Kite Mass Bottleneck"):
  - Rotor mass total : 11 kg (3 blades × 3.667 kg each)
  - TRPT shaft mass  : ~6.6 kg  (14 rings × 0.4 kg + 5 lines × 30 m × 0.006 kg/m)
  - Total airborne   : ~17.6 kg ≈ 17 kg

Tether geometry (DRR §5.2):
  - 5 Dyneema lines, 3 mm diameter, 30 m total length
  - Torque rings every 2 m → 14 intermediate rings  ((30/2) − 1 = 14)
  - Each ring ≈ 400 g (12 mm CFRP tubes + aluminium clevis connectors)

TRPT taper geometry (DRR Grasshopper FEA, rL Geometry Ratio: 0.740741):
  - trpt_hub_radius = 2.0 m  (ring radius at rotor end; self-consistent with rL=0.74,
    tether=30m, n_seg=15 → average L_seg=2.0m, r_bottom=0.96m)
  - trpt_rL_ratio = 0.74    (r/L geometry constraint per segment)

Aerodynamics (Rotor_TRTP_Sizing_Iteration2.xlsx AeroDyn BEM, NACA4412 profile):
  - Cp ≈ 0.22 at optimal TSR ≈ 4.1–4.2 (3-blade, 20°–30° elevation)
    [4kW: Cp=0.232, 7kW: Cp=0.223, 12kW: Cp=0.227 — consistent across sizes]
  - Previous Framework PDF value (0.15) was a conservative proxy; AeroDyn gives ≈0.22

Ground inertia (Mass Scaling PDF §"Drivetrain Mass and Inertia Matching"):
  - I_wheel = 0.019 kg·m², I_gen = 0.040 kg·m² → I_pto = 0.059 kg·m²
"""
function params_10kw()::SystemParams
    m_blade = 11.0 / 3.0   # 3 blades totalling 11 kg per Mass Scaling PDF
    m_ring  = 0.4           # ~400 g per ring (DRR §5.2: CFRP tubes + clevis connectors)
    i_pto   = 0.019 + 0.040 # wheel + generator (Mass Scaling PDF §"Drivetrain Mass")

    # Hub altitude = tether_length × sin(elevation_angle) = 30 × sin(30°) = 15 m
    # h_ref = hub altitude; v_wind_ref is rated hub wind speed (11 m/s, from DRR / AeroDyn data)
    return SystemParams(
        1.225,           # rho (kg/m³)
        11.0,            # v_wind_ref (m/s) — rated wind speed at h_ref (DRR; AeroDyn sizing)
        15.0,            # h_ref (m) — reference (measurement) altitude; set equal to hub altitude
                         #   (30 × sin(30°) = 15 m) so v_wind_ref is specified directly at hub.
                         #   Hellmann wind shear is active when elevation angle changes (hub
                         #   altitude moves away from h_ref). Use a met-mast height here if
                         #   you have wind data from a fixed anemometer at a different height.
        π / 6,           # elevation_angle = 30° (as physically built; DRR)
        deg2rad(70.0),   # lifter_elevation = 70° (typical operating angle for launch/landing lifter kite line)
        5.0,             # rotor_radius R (m) — Framework PDF §5.3 (aerodynamic outer radius)
        30.0,            # tether_length L₀ (m) — DRR §5.2 "For a 30m TRPT"
        2.0,             # trpt_hub_radius (m) — ring radius at rotor end; DRR Grasshopper rL=0.74,
                         #   tether=30m, n_seg=15 → avg L_seg=2.0m, r_bottom=0.96m
        0.74,            # trpt_rL_ratio — DRR Grasshopper "rL Geometry Ratio: 0.740741"
        5,               # n_lines — DRR §5.2 "5 tethers along the length"
        0.003,           # tether_diameter (m) — DRR §5.2: 3 mm Dyneema type 01505
        100e9,           # e_modulus (Pa) — Dyneema ~100 GPa
        14,              # n_rings — DRR §5.2: rings every 2 m → (30/2)−1 = 14
        m_ring,          # m_ring (kg)
        5,               # n_blades — one blade per tether line / polygon vertex (n_lines = 5)
        m_blade,         # m_blade (kg)
        0.22,            # cp — AeroDyn BEM (Rotor_TRTP_Sizing_Iteration2.xlsx); peak reference value
        i_pto,           # i_pto (kg·m²)
        11.0,            # k_mppt (N·m·s²/rad²) — MPPT gain: τ_gen = k × ω²
                         #   Derivation: ω_opt = λ_opt × v_rated / R = 4.1 × 11 / 5 = 9.02 rad/s
                         #   τ_net = τ_aero − τ_drag ≈ 889 N·m at rated operating point
                         #   k = τ_net / ω_opt² = 889 / 81.4 ≈ 10.9 → 11.0
                         #   Quadratic load law eliminates the bistability of the old linear c_pto
                         #   and gives correct MPPT at all wind speeds, not just rated.
        10_000.0,        # p_rated_w (W) — 10 kW rated
        deg2rad(23.0),   # β_min — blade-tip ground clearance floor
        deg2rad(67.0),   # β_max — lifter kite ceiling
        deg2rad(1.0),    # β_rate_max (rad/s) — 1 °/s actuation rate limit
        5e-5,            # kp_elev (rad/W/s) — 0.5 °/s per 1 kW overpower
    )
end

"""
    params_50kw()

50 kW target configuration, geometrically scaled from the 10 kW prototype.
Uses mass_scale() with the empirical 1.35 exponent from Mass Scaling PDF.
"""
function params_50kw()::SystemParams
    return mass_scale(params_10kw(), 10.0, 50.0)
end

"""
    mass_scale(base, base_power_kw, target_power_kw)

Scale a SystemParams to a new rated power using:
- Aerodynamic length scaling: x = (target/base)^(1/2) for all lengths.
  Rationale: P_aero ∝ R² (swept area), so R ∝ P^(1/2). Confirmed by AeroDyn BEM data:
  4 kW R=2.8 m → 12 kW R=4.8 m ≈ (12/4)^(1/2) × 2.8 ✓
- Empirical mass exponent 1.35 (Mass Scaling PDF §"The Empirical Mass Exponent"):
  m_scaled = m_base × (target/base)^1.35
- PTO inertia: I ∝ m × R² → exponent = 1.35 + 2×(1/2) = 2.35
- MPPT gain: k_mppt = τ/ω² where τ ∝ P^(3/2), ω² ∝ P^(-1) → k_mppt ∝ P^(5/2)
"""
function mass_scale(base::SystemParams,
                    base_power_kw::Float64,
                    target_power_kw::Float64)::SystemParams
    power_ratio = target_power_kw / base_power_kw
    geom_scale  = power_ratio^(1.0/2.0)   # linear dimension scale: R ∝ P^(1/2) since P ∝ R²
    mass_factor = power_ratio^1.35         # empirical mass exponent (Mass Scaling PDF)

    return SystemParams(
        base.rho,
        base.v_wind_ref,
        base.h_ref             * geom_scale,
        base.elevation_angle,                  # angle does not scale
        base.lifter_elevation,                 # angle does not scale
        base.rotor_radius      * geom_scale,
        base.tether_length     * geom_scale,
        base.trpt_hub_radius   * geom_scale,   # scales geometrically
        base.trpt_rL_ratio,                    # dimensionless ratio, does not scale
        base.n_lines,                          # topology does not scale
        base.tether_diameter   * geom_scale,
        base.e_modulus,                        # material property, unchanged
        base.n_rings,                          # topology does not scale
        base.m_ring            * mass_factor,
        base.n_blades,                         # topology does not scale
        base.m_blade           * mass_factor,
        base.cp,                               # aerodynamic constant, unchanged
        base.i_pto             * mass_factor * geom_scale^2,  # I ∝ m·R²: P^1.35 × P = P^2.35
        base.k_mppt            * power_ratio^2.5,             # k = τ/ω², τ ∝ P^(3/2), ω² ∝ P^(-1) → k ∝ P^(5/2)
        base.p_rated_w  * power_ratio,        # rated power scales linearly
        base.β_min,                            # angle does not scale
        base.β_max,                            # angle does not scale
        base.β_rate_max,                       # rate does not scale
        base.kp_elev    / power_ratio,         # larger turbine: same °/s per fractional overpower
    )
end
