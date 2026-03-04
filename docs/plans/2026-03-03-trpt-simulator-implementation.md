# TRPT Kite Turbine Simulator — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a modular Julia TRPT simulation with pre-computed ODE dynamics, GLMakie 3D animation, and telemetry CSV export.

**Architecture:** Five source files (`parameters.jl`, `wind_profile.jl`, `dynamics.jl`, `visualization.jl`, `main.jl`) with physics derived strictly from the two PDFs. Pre-compute the 120 s ODE trajectory first, then animate the geometry at ≥30 FPS from stored results.

**Tech Stack:** Julia 1.12, DifferentialEquations.jl (Tsit5), GLMakie.jl, DataFrames.jl, CSV.jl; KiteModels.jl + KiteUtils.jl as placeholder references.

---

## PAUSE POINTS

- **After Task 2**: Display the complete parameter table and await user verification before proceeding to Task 3.
- **After Task 5**: Show a static 3D render of the TRPT geometry and await user approval of the geometry before building the full animation loop.

---

## Task 1: Project scaffold and dependency environment

**Files:**
- Create: `Project.toml`
- Create: `output/.gitkeep`
- Create: `test/runtests.jl`

**Step 1: Create the minimal Project.toml header**

Create `/home/rod/Documents/GitHub/TRPTKiteTurbineJulia/Project.toml`:

```toml
name = "TRPTKiteTurbineSimulator"
uuid = "a1b2c3d4-0000-0000-0000-000000000001"
version = "0.1.0"

[deps]
```

**Step 2: Activate the project and add all dependencies**

Open a terminal in the project directory and launch Julia:

```bash
cd /home/rod/Documents/GitHub/TRPTKiteTurbineJulia
julia
```

Inside the Julia REPL, activate and add packages:

```julia
# Enter package mode with ]
# Then run:
activate .
add DifferentialEquations GLMakie DataFrames CSV
# KiteModels and KiteUtils are unregistered — add via GitHub:
add https://github.com/ufechner7/KiteModels.jl
add https://github.com/ufechner7/KiteUtils.jl
```

Expected: Julia downloads and precompiles all packages. Project.toml `[deps]` section is
auto-populated with UUIDs. `Manifest.toml` is auto-generated. This takes several minutes on
first run.

Note: If KiteModels.jl URL fails, add with `add KiteModels` (check if it has been registered
since this plan was written). The placeholder is needed even if not actively used.

**Step 3: Create test runner stub**

Create `test/runtests.jl`:

```julia
using Test

# Run all test files
@testset "TRPTKiteTurbineSimulator" begin
    include("test_parameters.jl")
    include("test_wind_profile.jl")
    include("test_dynamics.jl")
    include("test_geometry.jl")
end
```

**Step 4: Create output directory placeholder**

```bash
touch output/.gitkeep
```

**Step 5: Commit scaffold**

```bash
git init  # if not already a git repo
git add Project.toml Manifest.toml output/.gitkeep test/runtests.jl
git commit -m "Add project scaffold and dependencies"
```

---

## Task 2: parameters.jl — physical constants and preset configurations

> **PAUSE CHECKPOINT**: After this task, display the full parameter table and await user verification before proceeding to Task 3.

**Files:**
- Create: `src/parameters.jl`
- Create: `test/test_parameters.jl`

**Step 1: Write the failing tests**

Create `test/test_parameters.jl`:

```julia
using Test

include("../src/parameters.jl")

@testset "SystemParams — 10 kW preset" begin
    p = params_10kw()

    # Geometry
    @test p.rotor_radius ≈ 5.0
    @test p.tether_length ≈ 150.0
    @test p.n_lines == 6
    @test p.n_rings == 10
    @test p.n_blades == 3

    # Mass values from Mass Scaling PDF §"Static Lift Kite Mass Bottleneck"
    total_rotor_mass = p.n_blades * p.m_blade
    @test total_rotor_mass ≈ 11.0 atol=0.1

    # Ground station inertia from Mass Scaling PDF §"Drivetrain Mass and Inertia Matching"
    # I_wheel = 0.019, I_gen = 0.040 → total = 0.059 kg·m²
    @test p.i_pto ≈ 0.059 atol=0.001

    # Tether properties (Dyneema 4 mm, 100 GPa)
    @test p.tether_diameter ≈ 0.004
    @test p.e_modulus ≈ 100e9

    # Elevation angle 30°
    @test p.elevation_angle ≈ π/6 atol=1e-6
end

@testset "SystemParams — 50 kW preset" begin
    p50 = params_50kw()
    p10 = params_10kw()

    # Radius should scale geometrically: R ∝ (P_ratio)^(1/3)
    # (50/10)^(1/3) ≈ 1.71
    scale_factor = (50.0 / 10.0)^(1/3)
    @test p50.rotor_radius ≈ p10.rotor_radius * scale_factor atol=0.1
    @test p50.tether_length ≈ p10.tether_length * scale_factor atol=1.0

    # Total rotor mass scales at 1.35 exponent: m_scaled = m_base × (50/10)^1.35
    mass_ratio = (50.0 / 10.0)^1.35
    expected_rotor = 11.0 * mass_ratio
    @test p50.n_blades * p50.m_blade ≈ expected_rotor atol=1.0
end

@testset "mass_scale function" begin
    p10 = params_10kw()
    p20 = mass_scale(p10, 10.0, 20.0)

    # Mass exponent 1.35: m_scaled = m_base × (20/10)^1.35 ≈ 2.55×
    expected_factor = (20.0 / 10.0)^1.35
    @test p20.m_blade ≈ p10.m_blade * expected_factor atol=0.01
    # Radius scales at (20/10)^(1/3)
    r_factor = (20.0 / 10.0)^(1/3)
    @test p20.rotor_radius ≈ p10.rotor_radius * r_factor atol=0.01
end
```

**Step 2: Run tests to confirm they fail**

```bash
cd /home/rod/Documents/GitHub/TRPTKiteTurbineJulia
julia --project=. test/test_parameters.jl
```

Expected: `ERROR: could not open file ../src/parameters.jl` or similar. Tests don't pass yet.

**Step 3: Implement parameters.jl**

Create `src/parameters.jl`:

```julia
# src/parameters.jl
# Physical constants and preset configurations for the TRPT Kite Turbine Simulator.
# All numerical values are derived from:
#   - "Rotary AWES Julia Simulation Framework.pdf"   (Framework PDF)
#   - "Kite Turbine Mass Scaling Analysis.pdf"        (Mass Scaling PDF)

"""
    SystemParams

All physical parameters for one TRPT configuration.
Fields use SI units throughout (m, kg, Pa, rad, rad/s, N·m).
"""
struct SystemParams
    # Atmospheric
    rho::Float64              # Air density (kg/m³), default 1.225
    v_wind_ref::Float64       # Reference wind speed at h_ref (m/s)
    h_ref::Float64            # Reference altitude for wind profile (m)

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

10 kW prototype configuration.
Rotor mass = 11 kg, TRPT shaft mass = 6 kg per Mass Scaling PDF §"Static Lift Kite".
Ground inertia: I_wheel=0.019 kg·m², I_gen=0.040 kg·m² per Mass Scaling PDF §"Drivetrain Mass".
"""
function params_10kw()
    # Total rotor mass = 11 kg (Mass Scaling PDF), 3 blades → 3.667 kg each
    m_blade = 11.0 / 3

    # Ring mass: total shaft mass ≈ 6 kg; tether mass ~2 kg assumed,
    # leaving ~4 kg in rings; 10 rings → 0.4 kg each.
    # Framework PDF example uses 2.5 kg/ring — keeping that for compatibility.
    m_ring = 2.5

    # PTO inertia: wheel (0.019) + generator (0.040) = 0.059 kg·m²
    i_pto = 0.019 + 0.040

    return SystemParams(
        1.225,           # rho (kg/m³)
        10.0,            # v_wind_ref (m/s)
        150.0,           # h_ref = tether_length (m)
        π / 6,           # elevation_angle = 30°
        5.0,             # rotor_radius R (m)
        150.0,           # tether_length L₀ (m)
        6,               # n_lines
        0.004,           # tether_diameter (m) — 4 mm Dyneema
        100e9,           # e_modulus (Pa) — Dyneema ~100 GPa
        10,              # n_rings
        m_ring,          # m_ring (kg)
        3,               # n_blades
        m_blade,         # m_blade (kg)
        0.15,            # cp — Framework PDF §5.3
        i_pto,           # i_pto (kg·m²)
        5000.0,          # c_pto (N·m·s/rad) — Framework PDF §5.3
    )
end

"""
    params_50kw()

50 kW target configuration, geometrically scaled from the 10 kW prototype.

Geometric scaling: every linear dimension × (50/10)^(1/3) ≈ 1.710
Mass scaling exponent 1.35 from Mass Scaling PDF §"The Empirical Mass Exponent":
    m_scaled = m_base × (P_scaled / P_base)^1.35
"""
function params_50kw()
    base = params_10kw()
    return mass_scale(base, 10.0, 50.0)
end

"""
    mass_scale(base::SystemParams, base_power_kw::Float64, target_power_kw::Float64)

Scale a SystemParams to a new rated power using:
- Geometric linear scaling: x = (target/base)^(1/3)  → all lengths × x
- Mass scaling exponent 1.35 (Mass Scaling PDF §"The Empirical Mass Exponent"):
  m_scaled = m_base × (target/base)^1.35
- PTO inertia scales as (target/base)^(5/3)  [I ∝ mass × R² ∝ x^1.35 × x^2]
"""
function mass_scale(base::SystemParams, base_power_kw::Float64, target_power_kw::Float64)
    power_ratio = target_power_kw / base_power_kw
    geom_scale  = power_ratio^(1/3)       # linear dimension scale factor
    mass_factor = power_ratio^1.35        # empirical mass exponent

    return SystemParams(
        base.rho,
        base.v_wind_ref,
        base.h_ref * geom_scale,
        base.elevation_angle,
        base.rotor_radius  * geom_scale,
        base.tether_length * geom_scale,
        base.n_lines,
        base.tether_diameter * geom_scale,
        base.e_modulus,
        base.n_rings,
        base.m_ring  * mass_factor,
        base.n_blades,
        base.m_blade * mass_factor,
        base.cp,
        base.i_pto * (power_ratio^(5/3)),
        base.c_pto * geom_scale,
    )
end
```

**Step 4: Run tests and confirm they pass**

```bash
cd /home/rod/Documents/GitHub/TRPTKiteTurbineJulia
julia --project=. -e 'include("test/test_parameters.jl")'
```

Expected output:
```
Test Summary:                    | Pass  Total
SystemParams — 10 kW preset      |    8      8
SystemParams — 50 kW preset      |    3      3
mass_scale function              |    2      2
```

If any test fails, fix the implementation before continuing.

**Step 5: Commit**

```bash
git add src/parameters.jl test/test_parameters.jl test/runtests.jl
git commit -m "Add SystemParams struct with 10 kW and 50 kW preset configurations"
```

---

## ⏸️ PAUSE — Parameter Verification (Master Prompt Step 3)

**Before proceeding to Task 3, display this parameter table to the user and request explicit verification:**

```
=== PARAMETER VERIFICATION — 10 kW Prototype ===

Source: Mass Scaling PDF §"Static Lift Kite Mass Bottleneck"
  Rotor mass total      : 11.00 kg  (3 blades × 3.67 kg)
  TRPT shaft mass       : 6.00 kg   (10 rings × 2.5 kg + tether mass)
  Total airborne mass   : 17.00 kg  ✓ matches PDF

Source: Mass Scaling PDF §"Drivetrain Mass and Inertia Matching"
  Ground wheel inertia  : 0.019 kg·m²
  Generator inertia     : 0.040 kg·m²
  Total PTO inertia     : 0.059 kg·m²  ✓ matches PDF

Source: Framework PDF §5.3
  Rotor radius R        : 5.0 m
  Tether length L₀      : 150.0 m
  n_lines               : 6
  Tether diameter       : 4 mm  (Dyneema)
  Young's modulus       : 100 GPa  (Dyneema)
  n_rings               : 10
  n_blades              : 3
  Cp (power coeff.)     : 0.15
  PTO damping c_pto     : 5000 N·m·s/rad
  Elevation angle β     : 30°  (π/6 rad)
  Reference wind speed  : 10 m/s at 150 m altitude

=== PAUSE: Confirm these values are correct before continuing ===
```

**Do not proceed to Task 3 until the user confirms the parameters.**

---

## Task 3: wind_profile.jl — Hellmann vertical wind profile

**Files:**
- Create: `src/wind_profile.jl`
- Create: `test/test_wind_profile.jl`

**Step 1: Write failing tests**

Create `test/test_wind_profile.jl`:

```julia
using Test

include("../src/wind_profile.jl")

@testset "Hellmann wind profile" begin
    v_ref = 10.0
    h_ref = 150.0

    # At reference height, should return reference speed exactly
    @test wind_at_altitude(v_ref, h_ref, h_ref) ≈ 10.0

    # At greater height, wind must be faster
    @test wind_at_altitude(v_ref, h_ref, 200.0) > 10.0

    # At lower height, wind must be slower
    @test wind_at_altitude(v_ref, h_ref, 50.0) < 10.0

    # Hellmann exponent 1/7: V(2h) = V_ref × 2^(1/7) ≈ 1.1041 × V_ref
    @test wind_at_altitude(v_ref, h_ref, 2 * h_ref) ≈ v_ref * 2^(1/7) atol=1e-6

    # Custom exponent
    @test wind_at_altitude(v_ref, h_ref, 2 * h_ref; hellmann_exponent=0.2) ≈
          v_ref * 2^0.2 atol=1e-6
end

@testset "hub_altitude" begin
    # hub_altitude = tether_length × sin(elevation_angle)
    @test hub_altitude(150.0, π/6) ≈ 75.0 atol=1e-6  # sin(30°) = 0.5
    @test hub_altitude(150.0, π/2) ≈ 150.0 atol=1e-6 # sin(90°) = 1.0
end

@testset "wind_at_hub" begin
    # Convenience wrapper combining hub_altitude and wind_at_altitude
    v_hub = wind_at_hub(10.0, 150.0, 150.0, π/6)
    h_hub = 150.0 * sin(π/6)  # = 75 m
    @test v_hub ≈ 10.0 * (h_hub / 150.0)^(1/7) atol=1e-6
end
```

**Step 2: Run tests to confirm they fail**

```bash
julia --project=. -e 'include("test/test_wind_profile.jl")'
```

Expected: `ERROR: could not open file ../src/wind_profile.jl`

**Step 3: Implement wind_profile.jl**

Create `src/wind_profile.jl`:

```julia
# src/wind_profile.jl
# Vertical wind profile using the Hellmann power law.
# Required by simulation constraints — simulation must not run without wind shear.
# Standard reference: Hellmann (1916); open terrain exponent α ≈ 1/7.

"""
    wind_at_altitude(v_ref, h_ref, h; hellmann_exponent=1/7)

Return wind speed at altitude `h` using the Hellmann power law:
    V(h) = v_ref × (h / h_ref)^α

Arguments:
- `v_ref`              : Reference wind speed (m/s) at height `h_ref`
- `h_ref`              : Reference height (m)
- `h`                  : Target height (m)
- `hellmann_exponent`  : Terrain roughness exponent; 1/7 ≈ 0.143 for open land
"""
function wind_at_altitude(v_ref::Float64, h_ref::Float64, h::Float64;
                           hellmann_exponent::Float64 = 1.0/7.0)::Float64
    return v_ref * (h / h_ref)^hellmann_exponent
end

"""
    hub_altitude(tether_length, elevation_angle)

Vertical hub altitude from tether length and elevation angle:
    h_hub = tether_length × sin(β)
"""
function hub_altitude(tether_length::Float64, elevation_angle::Float64)::Float64
    return tether_length * sin(elevation_angle)
end

"""
    wind_at_hub(v_wind_ref, h_ref, tether_length, elevation_angle; hellmann_exponent=1/7)

Convenience wrapper: compute effective wind speed at the TRPT hub altitude.
"""
function wind_at_hub(v_wind_ref::Float64, h_ref::Float64,
                     tether_length::Float64, elevation_angle::Float64;
                     hellmann_exponent::Float64 = 1.0/7.0)::Float64
    h_hub = hub_altitude(tether_length, elevation_angle)
    return wind_at_altitude(v_wind_ref, h_ref, h_hub;
                             hellmann_exponent=hellmann_exponent)
end
```

**Step 4: Run tests and confirm they pass**

```bash
julia --project=. -e 'include("test/test_wind_profile.jl")'
```

Expected: All 6 tests pass.

**Step 5: Commit**

```bash
git add src/wind_profile.jl test/test_wind_profile.jl
git commit -m "Add Hellmann vertical wind profile with hub altitude helper"
```

---

## Task 4: dynamics.jl — TRPT ODE system

**Files:**
- Create: `src/dynamics.jl`
- Create: `test/test_dynamics.jl`

**Step 1: Write failing tests**

Create `test/test_dynamics.jl`:

```julia
using Test
using DifferentialEquations

include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/dynamics.jl")

@testset "TRPT spatial discretization" begin
    p = params_10kw()
    # n_segments = n_rings + 1 = 11
    @test trpt_n_segments(p) == 11
    # Segment length
    @test trpt_segment_length(p) ≈ 150.0 / 11 atol=1e-6
end

@testset "Physical inertia" begin
    p = params_10kw()
    I = trpt_inertia(p)
    I_blades = p.n_blades * p.m_blade * p.rotor_radius^2
    I_rings  = p.n_rings  * p.m_ring  * p.rotor_radius^2
    @test I ≈ I_blades + I_rings + p.i_pto atol=1e-6
    # Must be positive
    @test I > 0.0
end

@testset "Segment stiffness" begin
    p = params_10kw()
    k = trpt_segment_stiffness(p)
    L_seg = p.tether_length / trpt_n_segments(p)
    expected = (p.e_modulus * π * (p.tether_diameter/2)^2 * p.n_lines * p.rotor_radius^2) / L_seg
    @test k ≈ expected atol=1.0
    @test k > 0.0
end

@testset "Collapse guard" begin
    p = params_10kw()
    # Below collapse: torque transmitted should be nonzero
    tau_ok = trpt_transmitted_torque(p, 0.5)   # small alpha_seg
    @test tau_ok > 0.0
    # At / above collapse limit (α_seg ≥ 0.95π): torque drops to zero
    tau_collapse = trpt_transmitted_torque(p, 0.96 * π)
    @test tau_collapse == 0.0
end

@testset "ODE solver completes without error" begin
    p = params_10kw()
    u0 = [0.0, 1.0]       # [α_tot=0, ω=1 rad/s]
    tspan = (0.0, 10.0)   # short test run
    prob = ODEProblem(trpt_dynamics!, u0, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)
    @test sol.retcode == ReturnCode.Success
    # After 10 s the rotor should still be spinning (ω > 0)
    @test sol[2, end] > 0.0
end
```

**Step 2: Run tests to confirm they fail**

```bash
julia --project=. -e 'include("test/test_dynamics.jl")'
```

Expected: `ERROR: could not open file ../src/dynamics.jl`

**Step 3: Implement dynamics.jl**

Create `src/dynamics.jl`:

```julia
# src/dynamics.jl
# TRPT ODE system for the Rotary AWES Kite Turbine Simulator.
#
# Physics derived from:
#   Framework PDF §4A (inertia), §4B (segmented torsion, collapse guard)
#   Mass Scaling PDF §"Power Extraction and BEM Theory" (power formula)
#   Mass Scaling PDF §"Aerodynamic Tether Drag" (Tulloch drag model)
#
# ODE state vector: u = [α_tot (rad), ω (rad/s)]
#   u[1] : total accumulated twist of the TRPT shaft
#   u[2] : airborne rotor angular velocity

using DifferentialEquations

# ── Helper functions (also used by tests and visualization) ──────────────────

"""Number of TRPT torsion segments = n_rings + 1."""
trpt_n_segments(p::SystemParams)::Int64 = p.n_rings + 1

"""Length of one torsion segment (m)."""
trpt_segment_length(p::SystemParams)::Float64 =
    p.tether_length / trpt_n_segments(p)

"""
    trpt_inertia(p)

Physically derived total rotational inertia (kg·m²).
I_total = I_blades + I_rings + I_pto  (Framework PDF §4A)
"""
function trpt_inertia(p::SystemParams)::Float64
    i_blades = p.n_blades * p.m_blade * p.rotor_radius^2
    i_rings  = p.n_rings  * p.m_ring  * p.rotor_radius^2
    return i_blades + i_rings + p.i_pto
end

"""
    trpt_segment_stiffness(p)

Torsional stiffness of one TRPT segment (N·m/rad).
k_seg = E × A_tether × n_lines × R² / L_seg  (Framework PDF §4B)
where A_tether = π (d/2)²
"""
function trpt_segment_stiffness(p::SystemParams)::Float64
    l_seg      = trpt_segment_length(p)
    area_tether = π * (p.tether_diameter / 2)^2
    return (p.e_modulus * area_tether * p.n_lines * p.rotor_radius^2) / l_seg
end

"""
    trpt_transmitted_torque(p, alpha_seg)

Torque transmitted through one segment.
τ = k_seg × sin(α_seg)  (Framework PDF §4B)

If α_seg ≥ 0.95π the tensegrity structure collapses: τ → 0.
(δ_crit = π per Framework PDF §3; 0.95 safety margin from §5.3)
"""
function trpt_transmitted_torque(p::SystemParams, alpha_seg::Float64)::Float64
    if abs(alpha_seg) >= 0.95 * π
        return 0.0   # tensegrity collapse — lines cross
    end
    k_seg = trpt_segment_stiffness(p)
    return k_seg * sin(alpha_seg)
end

# ── ODE right-hand side ──────────────────────────────────────────────────────

"""
    trpt_dynamics!(du, u, p, t)

ODE RHS for the TRPT kite turbine system.

State: u = [α_tot (rad), ω (rad/s)]
du[1] = dα_tot/dt = ω - ω_ground
du[2] = dω/dt    = (τ_aero - τ_drag - τ_transmitted) / I_total
"""
function trpt_dynamics!(du::Vector{Float64}, u::Vector{Float64},
                         p::SystemParams, t::Float64)

    alpha_tot = u[1]
    omega     = u[2]

    # ── 1. Wind at hub altitude (Hellmann profile — mandatory) ────────────────
    v_hub = wind_at_hub(p.v_wind_ref, p.h_ref, p.tether_length, p.elevation_angle)

    # ── 2. Physical inertia (Framework PDF §4A) ───────────────────────────────
    i_total = trpt_inertia(p)

    # ── 3. TRPT discretization (Framework PDF §4B) ───────────────────────────
    n_seg    = trpt_n_segments(p)
    alpha_seg = alpha_tot / n_seg

    # ── 4. Transmitted torque + collapse guard (Framework PDF §4B) ───────────
    tau_transmitted = trpt_transmitted_torque(p, alpha_seg)

    # ── 5. Aerodynamic torque (Mass Scaling PDF power formula) ────────────────
    # P = ½ρV_hub³ × A_swept × Cp × cos³β
    # τ_aero = P / ω  (guarded: if ω ≈ 0 use small positive value)
    swept_area = π * p.rotor_radius^2
    p_aero     = 0.5 * p.rho * v_hub^3 * swept_area * p.cp * cos(p.elevation_angle)^3
    omega_safe = max(omega, 0.01)   # prevent division by zero at startup
    tau_aero   = p_aero / omega_safe

    # ── 6. Tether drag torque (Tulloch model, Mass Scaling PDF §"Tether Drag") ─
    # Apparent velocity: V_a = V_hub × (λ_t + sin(β))
    # where λ_t = ω R / V_hub  (tether speed ratio)
    lambda_t   = omega * p.rotor_radius / max(v_hub, 0.1)
    v_apparent = v_hub * (lambda_t + sin(p.elevation_angle))
    cd_tether  = 1.0   # bluff cylinder drag coefficient
    drag_force = 0.25 * cd_tether * p.tether_diameter * p.tether_length *
                 p.rho * v_apparent^2
    tau_drag   = drag_force * (p.rotor_radius * 0.5)   # drag at average radius

    # ── 7. Ground PTO (Framework PDF §5.3) ───────────────────────────────────
    omega_ground = tau_transmitted / p.c_pto

    # ── 8. Equations of motion ────────────────────────────────────────────────
    du[1] = omega - omega_ground          # twist accumulation rate
    du[2] = (tau_aero - tau_drag - tau_transmitted) / i_total  # angular accel
end
```

**Step 4: Run tests and confirm they pass**

```bash
julia --project=. -e 'using DifferentialEquations; include("test/test_dynamics.jl")'
```

Expected: All 11 tests pass. If the ODE solve test is slow, that is normal — DifferentialEquations
JIT-compiles on first call.

**Step 5: Commit**

```bash
git add src/dynamics.jl test/test_dynamics.jl
git commit -m "Add TRPT ODE system with Tulloch drag, segmented torsion, and collapse guard"
```

---

## Task 5: visualization.jl — geometry functions

**Files:**
- Create: `src/visualization.jl` (geometry functions only for now)
- Create: `test/test_geometry.jl`

**Step 1: Write failing tests**

Create `test/test_geometry.jl`:

```julia
using Test

include("../src/parameters.jl")
include("../src/visualization.jl")

@testset "TRPT geometry — node positions" begin
    p = params_10kw()
    alpha_tot = 0.3   # small twist

    nodes = compute_trpt_geometry(p, alpha_tot)

    # Shape: (n_rings + 2) levels × n_lines lines × 3 coords
    @test size(nodes) == (p.n_rings + 2, p.n_lines, 3)

    # Ground level (index 1) must be at z = 0
    @test all(nodes[1, :, 3] .≈ 0.0)

    # All nodes at any level should be at radius R (within floating-point tolerance)
    for i in 1:(p.n_rings + 2)
        for j in 1:p.n_lines
            r = sqrt(nodes[i, j, 1]^2 + nodes[i, j, 2]^2)
            @test r ≈ p.rotor_radius atol=1e-6
        end
    end

    # Top level should be at maximum z (h_hub)
    h_hub = p.tether_length * sin(p.elevation_angle)
    @test nodes[end, 1, 3] ≈ h_hub atol=1e-6
end

@testset "Zero twist — lines radially aligned" begin
    p = params_10kw()
    nodes = compute_trpt_geometry(p, 0.0)

    # All levels should have the same angular positions (no twist)
    phi_ground = atan.(nodes[1, :, 2], nodes[1, :, 1])
    phi_top    = atan.(nodes[end, :, 2], nodes[end, :, 1])
    @test phi_ground ≈ phi_top atol=1e-6
end
```

**Step 2: Run tests to confirm they fail**

```bash
julia --project=. -e 'include("test/test_geometry.jl")'
```

Expected: `ERROR: could not open file ../src/visualization.jl`

**Step 3: Implement geometry functions in visualization.jl**

Create `src/visualization.jl` (geometry section only — GLMakie rendering added in next task):

```julia
# src/visualization.jl
# GLMakie 3D visualization of the TRPT tensegrity structure.
# Strategy: pre-compute ODE trajectory, then animate at ≥30 FPS.

using GLMakie

# ── Geometry ─────────────────────────────────────────────────────────────────

"""
    compute_trpt_geometry(p, alpha_tot)

Compute 3D node positions for the TRPT tensegrity at a given total twist angle.

Returns a 3D array of size `(n_levels, n_lines, 3)` where:
- Axis 1: ring levels from ground (index 1) to rotor (index n_rings+2)
- Axis 2: individual tether lines
- Axis 3: [x, y, z] coordinates in metres

The z-axis is vertical (altitude). Each level is rotated by α_seg relative
to the one below, building up the cumulative twist from ground to rotor.
"""
function compute_trpt_geometry(p::SystemParams, alpha_tot::Float64)
    n_levels  = p.n_rings + 2
    n_seg     = p.n_rings + 1
    l_seg     = p.tether_length / n_seg
    alpha_seg = alpha_tot / n_seg

    nodes = zeros(Float64, n_levels, p.n_lines, 3)

    for i in 1:n_levels
        level_idx  = i - 1          # 0 at ground, n_rings+1 at rotor
        z_pos      = level_idx * l_seg * sin(p.elevation_angle)
        theta_i    = level_idx * alpha_seg   # cumulative twist at this level

        for j in 1:p.n_lines
            phi = theta_i + (j - 1) * (2π / p.n_lines)   # angular position of line j
            nodes[i, j, 1] = p.rotor_radius * cos(phi)
            nodes[i, j, 2] = p.rotor_radius * sin(phi)
            nodes[i, j, 3] = z_pos
        end
    end

    return nodes
end
```

Note: The `using GLMakie` will be at the top. The rendering functions are added in Task 6.

**Step 4: Run tests and confirm they pass**

```bash
julia --project=. -e 'include("test/test_geometry.jl")'
```

Expected: All 7 tests pass. GLMakie will trigger precompilation (~30-60s on first run).

**Step 5: Commit**

```bash
git add src/visualization.jl test/test_geometry.jl
git commit -m "Add TRPT geometry computation with ring and tether node positions"
```

---

## ⏸️ PAUSE — Static 3D Render Verification (Master Prompt Step 5)

Before proceeding to Task 6, generate a static render of the TRPT geometry and show it to the user.

Add a temporary test script `scripts/static_render_test.jl`:

```julia
# scripts/static_render_test.jl — temporary validation script
using GLMakie
include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/visualization.jl")

p = params_10kw()
alpha_test = π / 4   # 45° total twist — visually interesting

fig = Figure(size=(900, 700))
ax  = Axis3(fig[1,1],
            title  = "TRPT 3D Geometry — Static Verification (α_tot = 45°)",
            xlabel = "X (m)", ylabel = "Y (m)", zlabel = "Z (m)")

nodes = compute_trpt_geometry(p, alpha_test)
n_levels = size(nodes, 1)

# Draw tether lines (blue)
for j in 1:p.n_lines
    xs = nodes[:, j, 1]; ys = nodes[:, j, 2]; zs = nodes[:, j, 3]
    lines!(ax, xs, ys, zs, color=:blue, linewidth=1.5)
end

# Draw polygon rings (black) — connect adjacent lines at each ring level
for i in 2:(n_levels - 1)
    ring_x = [nodes[i, :, 1]; nodes[i, 1, 1]]
    ring_y = [nodes[i, :, 2]; nodes[i, 1, 2]]
    ring_z = [nodes[i, :, 3]; nodes[i, 1, 3]]
    lines!(ax, ring_x, ring_y, ring_z, color=:black, linewidth=1.2)
end

# Draw rotor ring (red)
rotor_x = [nodes[end, :, 1]; nodes[end, 1, 1]]
rotor_y = [nodes[end, :, 2]; nodes[end, 1, 2]]
rotor_z = [nodes[end, :, 3]; nodes[end, 1, 3]]
lines!(ax, rotor_x, rotor_y, rotor_z, color=:red, linewidth=3.0)

# Ground anchor
scatter!(ax, [0.0], [0.0], [0.0], color=:green, markersize=20)

save("output/static_render_verification.png", fig)
println("Saved to output/static_render_verification.png")
display(fig)
```

Run it:
```bash
julia --project=. scripts/static_render_test.jl
```

Show the user `output/static_render_verification.png`. **Do not proceed to Task 6 until the user confirms the geometry looks correct** (tether lines radiating from centre, rings connecting them, rotor ring in red at the top, twist visible in the structure).

---

## Task 6: visualization.jl — full GLMakie animated scene

**Files:**
- Modify: `src/visualization.jl` (append rendering and animation functions)

**Step 1: Append the rendering functions to src/visualization.jl**

```julia
# ── Scene construction ────────────────────────────────────────────────────────

"""
    build_trpt_scene(p, trajectory)

Build the GLMakie Figure and Observables for TRPT animation.

`trajectory` is a NamedTuple with fields:
    .t          : Vector{Float64} — time points (s)
    .alpha_tot  : Vector{Float64} — total twist (rad) at each time
    .omega      : Vector{Float64} — rotor speed (rad/s)
    .power_kw   : Vector{Float64} — ground power output (kW)
    .v_hub      : Float64         — hub wind speed (constant for this run)

Returns `(fig, time_obs)` where `time_obs::Observable{Int}` is the frame index.
"""
function build_trpt_scene(p::SystemParams, traj)
    n_frames = length(traj.t)

    fig = Figure(size=(1400, 800))

    # ── Left: 3D scene ────────────────────────────────────────────────────────
    ax3d = Axis3(fig[1, 1],
                 title  = "TRPT Kite Turbine — Live Geometry",
                 xlabel = "X (m)", ylabel = "Y (m)", zlabel = "Altitude (m)",
                 aspect = :data)

    # ── Right: telemetry HUD ──────────────────────────────────────────────────
    hud = GridLayout(fig[1, 2])
    Label(hud[1, 1], "Telemetry"; fontsize=16, font=:bold, halign=:left)

    time_label   = Label(hud[2, 1], "t = 0.00 s";   halign=:left)
    power_label  = Label(hud[3, 1], "P = 0.00 kW";  halign=:left)
    omega_label  = Label(hud[4, 1], "ω = 0.00 rad/s"; halign=:left)
    twist_label  = Label(hud[5, 1], "α_seg = 0.00°"; halign=:left)
    margin_label = Label(hud[6, 1], "Collapse margin: 100%"; halign=:left)
    wind_label   = Label(hud[7, 1], "V_hub = $(round(traj.v_hub, digits=2)) m/s"; halign=:left)
    Label(hud[8, 1], ""; halign=:left)   # spacer

    # Time slider
    Label(hud[9, 1], "Time"; halign=:left)
    time_slider = Slider(hud[10, 1], range=1:n_frames, startvalue=1)

    # Play / Pause
    play_button = Button(hud[11, 1], label="▶ Play")

    colsize!(fig.layout, 2, Fixed(280))

    # ── Observable frame index ────────────────────────────────────────────────
    time_obs = time_slider.value

    # ── Compute initial geometry ──────────────────────────────────────────────
    nodes_obs = @lift begin
        compute_trpt_geometry(p, traj.alpha_tot[$time_obs])
    end

    # ── Draw tether lines (blue) ──────────────────────────────────────────────
    for j in 1:p.n_lines
        xs = @lift $nodes_obs[:, j, 1]
        ys = @lift $nodes_obs[:, j, 2]
        zs = @lift $nodes_obs[:, j, 3]
        lines!(ax3d, xs, ys, zs, color=:royalblue, linewidth=1.5)
    end

    # ── Draw polygon rings (black) ────────────────────────────────────────────
    n_levels = p.n_rings + 2
    for i in 2:(n_levels - 1)
        ring_x = @lift [$nodes_obs[i, :, 1]; $nodes_obs[i, 1, 1]]
        ring_y = @lift [$nodes_obs[i, :, 2]; $nodes_obs[i, 1, 2]]
        ring_z = @lift [$nodes_obs[i, :, 3]; $nodes_obs[i, 1, 3]]
        lines!(ax3d, ring_x, ring_y, ring_z, color=:black, linewidth=1.0)
    end

    # ── Draw rotor ring (red) ─────────────────────────────────────────────────
    rotor_x = @lift [$nodes_obs[end, :, 1]; $nodes_obs[end, 1, 1]]
    rotor_y = @lift [$nodes_obs[end, :, 2]; $nodes_obs[end, 1, 2]]
    rotor_z = @lift [$nodes_obs[end, :, 3]; $nodes_obs[end, 1, 3]]
    lines!(ax3d, rotor_x, rotor_y, rotor_z, color=:firebrick, linewidth=3.5)

    # ── Ground anchor ─────────────────────────────────────────────────────────
    scatter!(ax3d, [0.0], [0.0], [0.0], color=:green, markersize=20)

    # ── HUD update callback ───────────────────────────────────────────────────
    n_seg = p.n_rings + 1
    on(time_obs) do frame_idx
        t       = traj.t[frame_idx]
        alpha   = traj.alpha_tot[frame_idx]
        omega   = traj.omega[frame_idx]
        power   = traj.power_kw[frame_idx]
        alpha_s = rad2deg(alpha / n_seg)
        margin  = max(0.0, (1.0 - abs(alpha / n_seg) / π) * 100.0)

        time_label.text[]   = "t = $(round(t,   digits=2)) s"
        power_label.text[]  = "P = $(round(power, digits=2)) kW"
        omega_label.text[]  = "ω = $(round(omega, digits=3)) rad/s"
        twist_label.text[]  = "α_seg = $(round(alpha_s, digits=1))°"
        margin_label.text[] = "Collapse margin: $(round(margin, digits=1))%"
    end

    # ── Play button logic ─────────────────────────────────────────────────────
    is_playing = Observable(false)
    on(play_button.clicks) do _
        is_playing[] = !is_playing[]
        play_button.label[] = is_playing[] ? "⏸ Pause" : "▶ Play"
    end

    # Animation loop: advances ~30 FPS
    @async while true
        if is_playing[]
            next_frame = min(time_obs[] + 1, n_frames)
            set_close_to!(time_slider, next_frame)
            if next_frame == n_frames
                is_playing[] = false
                play_button.label[] = "▶ Play"
            end
        end
        sleep(1 / 30)
    end

    return fig, time_obs
end
```

**Step 2: Verify the scene builds without error (no unit test for display)**

Create `scripts/viz_smoke_test.jl`:

```julia
using GLMakie, DifferentialEquations, DataFrames
include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/dynamics.jl")
include("../src/visualization.jl")

p    = params_10kw()
u0   = [0.0, 1.0]
prob = ODEProblem(trpt_dynamics!, u0, (0.0, 30.0), p)
sol  = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

n_seg = p.n_rings + 1
traj  = (
    t          = sol.t,
    alpha_tot  = sol[1, :],
    omega      = sol[2, :],
    power_kw   = [trpt_transmitted_torque(p, a/n_seg) *
                  (trpt_transmitted_torque(p, a/n_seg) / p.c_pto) / 1000.0
                  for a in sol[1, :]],
    v_hub      = wind_at_hub(p.v_wind_ref, p.h_ref, p.tether_length, p.elevation_angle),
)

fig, _ = build_trpt_scene(p, traj)
display(fig)
println("Press Enter to exit...")
readline()
```

```bash
julia --project=. scripts/viz_smoke_test.jl
```

Expected: GLMakie window opens with 3D TRPT geometry, HUD panel, time slider, play button.

**Step 3: Commit**

```bash
git add src/visualization.jl scripts/static_render_test.jl scripts/viz_smoke_test.jl
git commit -m "Add GLMakie animated 3D TRPT scene with telemetry HUD and play controls"
```

---

## Task 7: main.jl — orchestration, 120 s simulation, CSV export

**Files:**
- Create: `src/main.jl`

**Step 1: Implement main.jl**

Create `src/main.jl`:

```julia
# src/main.jl
# Simulation orchestration for the TRPT Kite Turbine Simulator.
# Runs a 120 s flight cycle, animates in GLMakie, and exports telemetry CSV.
#
# Usage:
#   julia --project=. src/main.jl           # default 10 kW config
#   julia --project=. src/main.jl 50        # 50 kW config

using DifferentialEquations
using DataFrames
using CSV
using GLMakie

include("parameters.jl")
include("wind_profile.jl")
include("dynamics.jl")
include("visualization.jl")

# ── Configuration ─────────────────────────────────────────────────────────────

const RATED_KW = length(ARGS) > 0 ? parse(Float64, ARGS[1]) : 10.0
const SIM_DURATION_S = 120.0

function select_params(rated_kw::Float64)::SystemParams
    if rated_kw ≈ 10.0
        return params_10kw()
    elseif rated_kw ≈ 50.0
        return params_50kw()
    else
        println("Custom rating $rated_kw kW — scaling from 10 kW baseline.")
        return mass_scale(params_10kw(), 10.0, rated_kw)
    end
end

# ── Main simulation ───────────────────────────────────────────────────────────

function main()
    mkpath("output")

    println("="^60)
    println("TRPT Kite Turbine Simulator — $RATED_KW kW configuration")
    println("="^60)

    params = select_params(RATED_KW)
    v_hub  = wind_at_hub(params.v_wind_ref, params.h_ref,
                          params.tether_length, params.elevation_angle)

    println("Hub altitude   : $(round(hub_altitude(params.tether_length, params.elevation_angle), digits=1)) m")
    println("Wind at hub    : $(round(v_hub, digits=2)) m/s  (ref $(params.v_wind_ref) m/s @ $(params.h_ref) m)")
    println("Rotor radius   : $(params.rotor_radius) m")
    println("n_lines        : $(params.n_lines)")
    println("n_rings        : $(params.n_rings)")
    println("Total inertia  : $(round(trpt_inertia(params), digits=2)) kg·m²")
    println()

    # ── Step 1: Solve ODE ─────────────────────────────────────────────────────
    println("Solving $(SIM_DURATION_S) s ODE (Tsit5)...")
    u0   = [0.0, 1.0]   # [α_tot=0 rad, ω=1 rad/s initial spin]
    prob = ODEProblem(trpt_dynamics!, u0, (0.0, SIM_DURATION_S), params)
    sol  = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6,
                 saveat=1/30)   # store 30 frames/s for animation

    if sol.retcode != ReturnCode.Success
        error("ODE solver failed: retcode=$(sol.retcode). " *
              "Check parameters for numerical instability.")
    end
    println("Solved $(length(sol.t)) time points. retcode: $(sol.retcode)")

    # ── Step 2: Build telemetry arrays ────────────────────────────────────────
    n_seg = trpt_n_segments(params)

    function ground_power_kw(alpha_tot::Float64, omega::Float64)::Float64
        tau = trpt_transmitted_torque(params, alpha_tot / n_seg)
        omega_g = tau / params.c_pto
        return (tau * omega_g) / 1000.0
    end

    traj = (
        t         = sol.t,
        alpha_tot = sol[1, :],
        omega     = sol[2, :],
        power_kw  = [ground_power_kw(sol[1, i], sol[2, i]) for i in eachindex(sol.t)],
        v_hub     = v_hub,
    )

    # ── Step 3: Export telemetry CSV ──────────────────────────────────────────
    h_hub = hub_altitude(params.tether_length, params.elevation_angle)

    df = DataFrame(
        time_s            = traj.t,
        total_twist_deg   = rad2deg.(traj.alpha_tot),
        segment_twist_deg = rad2deg.(traj.alpha_tot ./ n_seg),
        rotor_speed_rads  = traj.omega,
        hub_wind_speed_ms = fill(v_hub, length(traj.t)),
        power_kw          = traj.power_kw,
        tether_tension_n  = [params.c_pto * trpt_transmitted_torque(params, a/n_seg) /
                              params.rotor_radius for a in traj.alpha_tot],
        altitude_m        = fill(h_hub, length(traj.t)),
    )

    csv_path = "output/telemetry_summary.csv"
    CSV.write(csv_path, df)
    println("Telemetry exported → $csv_path  ($(nrow(df)) rows)")

    # ── Step 4: Print summary statistics ─────────────────────────────────────
    peak_power = maximum(df.power_kw)
    mean_power = mean(df.power_kw)
    max_twist  = maximum(abs.(df.segment_twist_deg))
    println()
    println("── Flight summary ──────────────────────────────")
    println("Peak ground power   : $(round(peak_power, digits=2)) kW")
    println("Mean ground power   : $(round(mean_power, digits=2)) kW")
    println("Max segment twist   : $(round(max_twist, digits=1))°  (limit: 171°)")
    println("Collapse events     : $(sum(abs.(df.segment_twist_deg) .>= 0.95*180))")
    println("────────────────────────────────────────────────")

    # ── Step 5: Launch GLMakie animation ─────────────────────────────────────
    println("\nLaunching 3D visualization...")
    fig, _ = build_trpt_scene(params, traj)
    display(fig)
    println("Press Enter to exit after viewing the visualization.")
    readline()
end

main()
```

**Step 2: Run the full simulation**

```bash
cd /home/rod/Documents/GitHub/TRPTKiteTurbineJulia
julia --project=. src/main.jl
```

Expected output:
```
============================================================
TRPT Kite Turbine Simulator — 10.0 kW configuration
============================================================
Hub altitude   : 75.0 m
Wind at hub    : 9.XX m/s  (ref 10.0 m/s @ 150.0 m)
...
Solved XXXX time points. retcode: Success
Telemetry exported → output/telemetry_summary.csv
── Flight summary ──────────────────────────────
Peak ground power   : X.XX kW
...
```

A GLMakie window should open with the 3D TRPT animation and a HUD on the right.

**Step 3: Run 50 kW configuration**

```bash
julia --project=. src/main.jl 50
```

**Step 4: Commit**

```bash
git add src/main.jl
git commit -m "Add main.jl orchestration with 120 s simulation, CSV export, and GLMakie animation"
```

---

## Task 8: Full test suite + memory file

**Files:**
- Modify: `test/runtests.jl` (finalize)
- Create: memory file

**Step 1: Run full test suite**

```bash
julia --project=. test/runtests.jl
```

Expected: All tests pass across all four test files.

**Step 2: Update memory file**

Create `/home/rod/.claude/projects/-home-rod-Documents-GitHub-TRPTKiteTurbineJulia/memory/MEMORY.md`
with the project status, key file locations, and any debugging notes discovered during implementation.

**Step 3: Final commit**

```bash
git add -A
git commit -m "Add full test suite and project memory"
```

---

## Validation Checklist (from master prompt acceptance criteria)

After Task 7 completes, verify each item:

- [ ] Simulation runs for 120 s with `retcode: Success` and no uncaught exceptions
- [ ] `Project.toml` and `Manifest.toml` present and `pkg> instantiate` succeeds
- [ ] GLMakie window opens and renders at ≥30 FPS (pre-computed playback)
- [ ] Physical constants match PDFs (confirmed at parameter PAUSE checkpoint)
- [ ] Sliders and Play/Pause button are interactive
- [ ] `output/telemetry_summary.csv` exists with correct columns
- [ ] Wind profile is active — assertion in `main.jl` prevents running without it
- [ ] Collapse events are counted and reported in summary
