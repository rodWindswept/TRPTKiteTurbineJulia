# TRPT Kite Turbine Simulator — Design Document

**Date:** 2026-03-03
**Author:** Rod Read / Claude Code
**Status:** Approved

---

## 1. Purpose

Build a modular Julia simulation environment for Rotary Airborne Wind Energy Systems (RAWES)
using Tensile Rotary Power Transmission (TRPT). The simulator solves the coupled rotational
dynamics of the airborne kite rotor and TRPT tensegrity shaft, produces a live 3D GLMakie
visualization, and exports a telemetry CSV for post-flight analysis.

Source-of-truth specifications:
- `Rotary AWES Julia Simulation Framework.pdf` — ODE structure, segmented torsion, parameters
- `Kite Turbine Mass Scaling Analysis.pdf` — empirical constants, tether drag model, mass scaling

---

## 2. Architecture

```
TRPTKiteTurbineJulia/
├── Project.toml
├── Manifest.toml
├── output/
│   └── .gitkeep
├── docs/plans/
│   └── 2026-03-03-trpt-simulator-design.md
└── src/
    ├── parameters.jl      # Physical constants, structs, preset configs, mass scaling
    ├── wind_profile.jl    # Hellmann power-law vertical wind profile
    ├── dynamics.jl        # TRPT ODE system, tether drag, segmented torsion
    ├── visualization.jl   # GLMakie 3D scene, animation, telemetry HUD
    └── main.jl            # Orchestration: solve → animate → export CSV
```

---

## 3. Module Specifications

### 3.1 parameters.jl

**Struct: `SystemParams`**

| Field | Type | Description | Source |
|---|---|---|---|
| `rho` | Float64 | Air density (kg/m³), default 1.225 | Framework PDF §3 |
| `v_wind_ref` | Float64 | Reference wind speed at h_ref (m/s) | Framework PDF §3 |
| `h_ref` | Float64 | Reference altitude for wind profile (m) | Wind profile |
| `elevation_angle` | Float64 | TRPT shaft elevation angle β (rad) | Mass Scaling PDF |
| `rotor_radius` | Float64 | Airborne ring radius R (m) | Framework PDF §3 |
| `n_lines` | Int64 | Number of TRPT tether lines | Framework PDF §3 |
| `tether_diameter` | Float64 | Tether line diameter (m) | Framework PDF §3 |
| `e_modulus` | Float64 | Tether Young's modulus (Pa), Dyneema ~100 GPa | Framework PDF §3 |
| `tether_length` | Float64 | Unstretched total TRPT length L₀ (m) | Framework PDF §3 |
| `n_rings` | Int64 | Number of polygon spacer rings | Framework PDF §3 |
| `m_ring` | Float64 | Mass per polygon ring (kg) | Framework PDF §3 |
| `n_blades` | Int64 | Number of lifting blades/kites | Framework PDF §3 |
| `m_blade` | Float64 | Mass per blade (kg) | Mass Scaling PDF |
| `cp` | Float64 | Rotor power coefficient, default 0.15 | Framework PDF §5.3 |
| `i_pto` | Float64 | Ground PTO inertia (kg·m²) | Mass Scaling PDF |
| `c_pto` | Float64 | PTO damping coefficient (N·m·s/rad) | Framework PDF §5.3 |

**Preset configurations:**

`params_10kw()` — 10 kW prototype, empirically validated:
- rotor_radius = 5.0 m, tether_length = 150 m, elevation_angle = π/6 (30°)
- n_lines = 6, tether_diameter = 0.004 m, e_modulus = 100e9 Pa
- n_rings = 10, m_ring = 2.5 kg (total shaft = 6 kg → 10 rings at some combination + tether mass)
- n_blades = 3, m_blade = 3.67 kg (total rotor = 11 kg per Mass Scaling PDF §"Static Lift Kite Mass Bottleneck")
- i_pto = 0.059 kg·m² (0.019 wheel + 0.040 generator, Mass Scaling PDF §"Drivetrain Mass")
- c_pto = 5000.0 N·m·s/rad, h_ref = 150.0 m

`params_50kw()` — 50 kW target, scaled from 10 kW:
- Scale factor x = (50/10)^(1/3) ≈ 1.71 (geometric scaling)
- rotor_radius = 5.0 × 1.71 ≈ 8.55 m
- tether_length = 150 × 1.71 ≈ 256 m
- m_blade scaled by mass exponent: m_scaled = m_base × (50/10)^1.35 ≈ 8.1× per unit
- i_pto scaled proportionally

**Function: `mass_scale(base::SystemParams, target_kw::Float64)`**
- Returns a new `SystemParams` with blade and ring masses scaled by
  `(target_kw / base_rated_kw)^1.35` per empirical exponent in Mass Scaling PDF.

---

### 3.2 wind_profile.jl

**Hellmann power law** (required by constraint — simulation must not run without it):

```
V(h) = V_ref × (h / h_ref)^α
```

where α = 1/7 ≈ 0.143 (standard Hellmann exponent for open land terrain).

Hub altitude is computed from tether length and elevation angle:
```
h_hub = tether_length × sin(elevation_angle)
```

**Function: `wind_at_altitude(v_ref, h_ref, h; hellmann_exponent=1/7)`**

Returns effective wind speed at the hub altitude. Called by `dynamics.jl` before every ODE
evaluation to ensure wind shear is always applied.

---

### 3.3 dynamics.jl

**ODE state vector:** `u = [α_tot, ω]`
- `u[1]` = total TRPT twist angle (rad)
- `u[2]` = airborne rotor angular velocity (rad/s)

**Physics pipeline (each timestep):**

1. **Wind at hub altitude** — via `wind_at_altitude()`, ensuring Hellmann profile.

2. **Physical inertia** (Framework PDF §4A):
   - `I_blades = n_blades × m_blade × R²`
   - `I_rings = n_rings × m_ring × R²`
   - `I_total = I_blades + I_rings + i_pto`

3. **TRPT spatial discretization** (Framework PDF §4B):
   - `n_segments = n_rings + 1`
   - `α_seg = α_tot / n_segments` ← twist between consecutive rings
   - `L_seg = tether_length / n_segments`

4. **Segmented torsional stiffness** (Framework PDF §4B):
   - `k_seg = (e_modulus × π × (d/2)² × n_lines × R²) / L_seg`
   - `τ_transmitted = k_seg × sin(α_seg)`

5. **Collapse guard** (Framework PDF §2, §4B):
   - If `|α_seg| ≥ 0.95π`: `τ_transmitted = 0` (tensegrity structural collapse)

6. **Aerodynamic torque** (Mass Scaling PDF power formula):
   - `P_aero = 0.5 × ρ × V_hub³ × π×R² × cp × cos³(β)`
   - `τ_aero = P_aero / ω` (guarded against ω → 0)

7. **Tether drag** (Tulloch model, Mass Scaling PDF §"Aerodynamic Tether Drag"):
   - Apparent velocity: `V_a = V_hub × (λ_t + sin(β))`
     where `λ_t = ω × R / V_hub` (tether speed ratio)
   - Averaged across the ring circumference (single representative phase angle)
   - `drag_force = 0.25 × cd_tether × tether_diameter × tether_length × ρ × V_a²`
   - `τ_drag = drag_force × R × 0.5`

8. **Ground PTO** (Framework PDF §5.3):
   - `ω_ground = τ_transmitted / c_pto`

9. **Equations of motion:**
   - `dα_tot/dt = ω - ω_ground`
   - `dω/dt = (τ_aero - τ_drag - τ_transmitted) / I_total`

**Solver:** `Tsit5()` with `reltol=1e-6, abstol=1e-6`. On solver instability the main script
logs the error and does not report "success" (per constraints).

---

### 3.4 visualization.jl

**Strategy:** Pre-compute the full 120 s ODE trajectory, then animate at ≥30 FPS.

**GLMakie scene layout:**
- Left panel (70% width): 3D TRPT geometry
  - Ground anchor: green sphere at origin
  - Tether lines: blue `linesegments!` connecting ring nodes
  - Polygon spacer rings: black lines connecting adjacent tether nodes at each level
  - Rotor ring + blade spokes: red thick lines at top
- Right panel (30% width): telemetry HUD + controls
  - Time slider (scrub through 120 s)
  - Play/Pause button
  - Live readouts: Power (kW), Segment Twist (°), ω (rad/s), Collapse margin (%)

**`compute_trpt_geometry(params, α_tot)`** → returns node coordinates matrix:
- For each ring level `i` in `0 : n_rings+1`:
  - z_pos = i × L_seg × sin(elevation_angle)
  - θ_i = i × α_seg
  - For each line `j`: `(R cos(θ_i + φ_j), R sin(θ_i + φ_j), z_pos)`
    where `φ_j = 2π(j-1)/n_lines`

**Observable connections:**
- `time_obs::Observable{Int}` drives geometry and HUD updates
- Slider updates `time_obs`; play loop increments `time_obs` at 30 FPS

---

### 3.5 main.jl

1. Load parameters (default: `params_10kw()`)
2. Validate wind profile is active (assertion)
3. Solve ODE: `solve(prob, Tsit5(), ...)`, check for solver errors
4. Animate via `visualization.jl`
5. Export `output/telemetry_summary.csv` with columns:
   - `time_s`, `total_twist_deg`, `segment_twist_deg`, `rotor_speed_rads`,
     `hub_wind_speed_ms`, `power_kw`, `tether_tension_n`, `altitude_m`

---

## 4. Dependencies (Project.toml)

| Package | Purpose |
|---|---|
| `DifferentialEquations` | ODE solver |
| `GLMakie` | Interactive 3D visualization |
| `DataFrames` | Tabular data handling |
| `CSV` | Telemetry export |
| `KiteModels` | Placeholder — respected physics reference; not used in TRPT physics yet |
| `KiteUtils` | Placeholder — may provide useful utilities in future |

---

## 5. Physical Constants from PDFs

| Constant | Value | Source |
|---|---|---|
| Rotor mass (10 kW) | 11 kg | Mass Scaling PDF §"Static Lift Kite Mass Bottleneck" |
| TRPT shaft mass (10 kW) | 6 kg | Mass Scaling PDF §"Static Lift Kite Mass Bottleneck" |
| Total airborne mass (10 kW) | 17 kg | Mass Scaling PDF §"Static Lift Kite Mass Bottleneck" |
| Ground wheel inertia | 0.019 kg·m² | Mass Scaling PDF §"Drivetrain Mass and Inertia Matching" |
| Generator inertia | 0.040 kg·m² | Mass Scaling PDF §"Drivetrain Mass and Inertia Matching" |
| Mass scaling exponent | 1.35 | Mass Scaling PDF §"The Empirical Mass Exponent" |
| Dyneema Young's modulus | 100 GPa | Framework PDF §5.3 |
| Default Cp | 0.15 | Framework PDF §5.3 |
| Hellmann exponent | 1/7 ≈ 0.143 | Standard (open terrain) |
| δ_crit per segment | π rad (180°) | Framework PDF §3, §4B |
| Cut-in wind speed | 4 m/s | Mass Scaling PDF §"Static Lift Kite Mass Bottleneck" |
| Safety wind limit | 26 m/s | Mass Scaling PDF §"Automated Control and the BackBot" |

---

## 6. Known Limitations (per Framework PDF §6)

1. **Rigid body assumption** — blade aeroelastic deformation not modelled
2. **Bulk Cp** — true BEM loop not implemented; Cp=0.15 is a reasonable proxy
3. **Lift kite decoupled** — TRPT axis assumed static; lift kite dynamics not modelled

These are documented limitations, not defects, and are consistent with the scope of
the Framework PDF.

---

## 7. Step PAUSE Points

Per the master prompt:
- **After `parameters.jl`**: User verifies the parameter table before dynamics are written
- **After `visualization.jl`**: User approves static 3D geometry before animation loop
