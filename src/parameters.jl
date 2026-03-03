# src/parameters.jl
# Physical constants and preset configurations for the TRPT Kite Turbine Simulator.
# All numerical values derived from:
#   - "Rotary AWES Julia Simulation Framework.pdf"  (Framework PDF)
#   - "Kite Turbine Mass Scaling Analysis.pdf"       (Mass Scaling PDF)

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
    rotor_radius::Float64     # Airborne ring radius R (m)
    tether_length::Float64    # Unstretched total TRPT length L₀ (m)

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
    cp::Float64               # Rotor power coefficient; Framework PDF §5.3 ≈ 0.15

    # Ground station — Mass Scaling PDF §"Drivetrain Mass and Inertia Matching"
    i_pto::Float64            # Total PTO rotational inertia (kg·m²)
    c_pto::Float64            # PTO damping coefficient (N·m·s/rad)
end

"""
    params_10kw()

10 kW prototype configuration using empirically validated values from the PDFs.

Mass budget (Mass Scaling PDF §"Static Lift Kite Mass Bottleneck"):
  - Rotor mass total : 11 kg (3 blades × 3.667 kg each)
  - TRPT shaft mass  : 6 kg
  - Total airborne   : 17 kg

Ground inertia (Mass Scaling PDF §"Drivetrain Mass and Inertia Matching"):
  - I_wheel = 0.019 kg·m², I_gen = 0.040 kg·m² → I_pto = 0.059 kg·m²
"""
function params_10kw()::SystemParams
    m_blade = 11.0 / 3.0   # 3 blades totalling 11 kg per Mass Scaling PDF
    m_ring  = 2.5           # Framework PDF §5.3 example value
    i_pto   = 0.019 + 0.040 # wheel + generator (Mass Scaling PDF §"Drivetrain Mass")

    return SystemParams(
        1.225,           # rho (kg/m³)
        10.0,            # v_wind_ref (m/s)
        150.0,           # h_ref (m) = tether length
        π / 6,           # elevation_angle = 30° per Mass Scaling PDF
        5.0,             # rotor_radius R (m) — Framework PDF §5.3
        150.0,           # tether_length L₀ (m)
        6,               # n_lines — Framework PDF §5.3
        0.004,           # tether_diameter (m) — 4 mm Dyneema
        100e9,           # e_modulus (Pa) — Dyneema ~100 GPa
        10,              # n_rings — Framework PDF §5.3
        m_ring,          # m_ring (kg)
        3,               # n_blades — Framework PDF §5.3
        m_blade,         # m_blade (kg)
        0.15,            # cp — Framework PDF §5.3
        i_pto,           # i_pto (kg·m²)
        5000.0,          # c_pto (N·m·s/rad) — Framework PDF §5.3
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
- Geometric linear scaling: x = (target/base)^(1/3) for all lengths
- Empirical mass exponent 1.35 (Mass Scaling PDF §"The Empirical Mass Exponent"):
  m_scaled = m_base × (target/base)^1.35
- PTO inertia: I ∝ mass × R² → exponent = 1.35 + 2/3 = 2.017 ≈ (1.35 + 2×1/3)
"""
function mass_scale(base::SystemParams,
                    base_power_kw::Float64,
                    target_power_kw::Float64)::SystemParams
    power_ratio = target_power_kw / base_power_kw
    geom_scale  = power_ratio^(1.0/3.0)   # linear dimension scale factor
    mass_factor = power_ratio^1.35         # empirical mass exponent

    return SystemParams(
        base.rho,
        base.v_wind_ref,
        base.h_ref             * geom_scale,
        base.elevation_angle,                  # angle does not scale
        base.rotor_radius      * geom_scale,
        base.tether_length     * geom_scale,
        base.n_lines,                          # topology does not scale
        base.tether_diameter   * geom_scale,
        base.e_modulus,                        # material property, unchanged
        base.n_rings,                          # topology does not scale
        base.m_ring            * mass_factor,
        base.n_blades,                         # topology does not scale
        base.m_blade           * mass_factor,
        base.cp,                               # aerodynamic constant, unchanged
        base.i_pto             * (power_ratio^(5.0/3.0)),  # I ∝ m·R²
        base.c_pto             * geom_scale,
    )
end
