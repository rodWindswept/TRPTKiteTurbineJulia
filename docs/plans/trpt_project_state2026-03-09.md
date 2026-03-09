# TRPTKiteTurbineJulia — State of Project & Development Priorities

*Analysis date: 2026-03-09*

---

## What Has Been Built

### Phase 1 — Core Simulator ✅ Complete

| Module | File | Status |
|--------|------|--------|
| Physical parameters + scaling | `src/parameters.jl` | ✅ 24-field SystemParams incl. limiter fields |
| Hellmann wind profile + wind models | `src/wind_profile.jl` | ✅ Complete |
| TSR-dependent aerodynamics | `src/aerodynamics.jl` | ✅ BEM Cp(λ)/CT(λ) lookup + interpolation |
| TRPT ODE dynamics | `src/dynamics.jl` | ✅ 2-state + 3-state (limiter) ODEs |
| Entry point / orchestration | `src/main.jl` | ✅ Complete |
| Test suite | `test/` (395 tests) | ✅ All pass |

### Phase 2 — High-Fidelity 3D Viz Engine ✅ Complete

| Module | File | Status |
|--------|------|--------|
| Inclined shaft 3D geometry + blades | `src/geometry.jl` | ✅ Implemented |
| Per-element force analysis (ForceState) | `src/force_analysis.jl` | ✅ Implemented |
| GLMakie scene + UI | `src/visualization.jl` | ✅ Implemented |

### Phase 3 — Elevation Power Limiter ✅ Complete

**Goal:** Proportional elevation-angle controller that raises β when P > P_rated, preventing runaway power above rated wind.

| Feature | Status |
|---------|--------|
| 5 limiter fields in SystemParams (p_rated_w, β_min, β_max, β_rate_max, kp_elev) | ✅ |
| 3-state ODEs: trpt_ode_limited!, trpt_ode_wind_limited! | ✅ |
| instantaneous_power/3 (3-arg variant for 3-state trajectories) | ✅ |
| β trace in traj_norm — traj_obs type always includes beta | ✅ |
| 3D shaft animates live β from trajectory during playback | ✅ |
| HUD Elevation β row — shows live β during scrub/play | ✅ |
| Re-run button fully working (was broken — missing 5 SystemParams fields) | ✅ fixed |
| Re-run warm-starts ODE at optimal TSR (was cold-start ω=1.0) | ✅ fixed |
| Elevation slider minimum clamped to 23° (β_min floor) | ✅ |
| run_wind_ramp_limited.jl — 6→23 m/s with active limiter | ✅ |
| power_curve.jl — 4 curves: MPPT + limited × 10 kW + 50 kW, sweep to 23 m/s | ✅ |
| run_wind_ramp.jl — extended to 23 m/s, warm-start | ✅ |

---

## Physics Parameters (10 kW prototype — confirmed from DRR §5.2 and AeroDyn BEM)

| Parameter | Value | Source |
|-----------|-------|--------|
| tether_length | 30 m | DRR §5.2 |
| n_lines | 5 | DRR §5.2 |
| tether_diameter | 3 mm (Dyneema 01505) | DRR §5.2 |
| n_rings | 14 (30/2 − 1) | DRR §5.2 |
| rotor_radius | 5.0 m | Framework PDF §5.3 |
| trpt_hub_radius | 2.0 m | DRR Grasshopper FEA |
| trpt_rL_ratio | 0.74 | DRR Grasshopper "rL = 0.740741" |
| cp | 0.22 (peak BEM) | AeroDyn NACA4412 3-blade |
| elevation_angle | 30° (π/6) | As physically built |
| k_mppt | 11.0 N·m·s²/rad² | τ_rated/ω_rated² = 889/9.02² ≈ 10.9 |
| p_rated_w | 10,000 W | Design target |
| β_min | 23° | Blade-tip ground clearance floor |
| β_max | 67° | Lifter kite pull limit |
| β_rate_max | 1°/s | Representative servo / lifter kite actuation rate |
| kp_elev | 5×10⁻⁵ rad/W/s | ≈ 0.5°/s per 1 kW overpower |

---

## Current Script Inventory

| Script | ODE | Wind | Output |
|--------|-----|------|--------|
| `run_dashboard_10kw.jl` | 2-state | steady 11 m/s | interactive dashboard |
| `run_dashboard_50kw.jl` | 2-state | steady | interactive dashboard |
| `run_wind_ramp.jl` | 2-state | 6→23 m/s | interactive dashboard (no limiter) |
| `run_wind_ramp_limited.jl` | **3-state** | 6→23 m/s | interactive dashboard (limiter active, β animates) |
| `run_gust.jl` | 2-state | 17 m/s Hann gust | interactive dashboard |
| `run_turbulent.jl` | 2-state | IEC Class A 15% TI | interactive dashboard |
| `record_simulation.jl` | 2-state | steady | MP4 video |
| `record_50kw.jl` | 2-state | steady | MP4 video |
| `power_curve.jl` | 2-state + **3-state** | 4–23 m/s sweep | CSV + PNG (4 curves) |

---

## Known Issues / Tech Debt

| Issue | Severity | Notes |
|-------|----------|-------|
| `@async` animation loop never exits | Low | Immortal task in interactive sessions; fine for script use |
| Duplicated r_top/r_bottom/k_eff block in trpt_ode! and instantaneous_power | Low | Candidate for extraction once stable |
| Force colour contrast low in dashboard | Medium | Pre-tension dominates; variable component is small relative to SWL |
| Bearing/lift kite line shown along shaft axis | Low/cosmetic | In real life lift kite is at 70° (p.lifter_elevation), separate from shaft — visualisation approximation |
| ⏸ pause glyph missing from Makie font on Linux | Fixed | Replaced with `\|\| Pause` (ASCII) |

---

## Physics Gaps (documented, not bugs)

1. **Rigid body** — no blade aeroelastic deformation
2. **Bulk Cp(λ) table** — BEM loop not run per-frame; Cp(λ) is interpolated from pre-computed tables
3. **Lift kite decoupled** — TRPT elevation angle β is a lumped control parameter; actual lift kite aerodynamics not modelled
4. **Blade coning / bank angle** — real blades are banked; swept disc is conical not flat; cos³(β) correction is for shaft inclination only, not coning
5. **Ring loading direction** — force_analysis.jl assumes compression; centrifugal loading at high ω may put rings in tension (different material sizing implications)

---

## Suggested Next Development Priorities

### Priority 1 — Lift kite geometry visualisation
*Low-risk, high-value visual improvement (plan already drafted — see previous session)*

Move the bearing to the rotor hub; show lift line leaving at `p.lifter_elevation = 70°` (not along shaft axis). Currently the bearing floats above the rotor along the shaft, misrepresenting the physical topology. The shaft (30°–67°) and lift kite line (70°) should form a visible angle in the scene.

### Priority 2 — Turbulent and gust scripts: add limited versions
*Extend existing scripts to use 3-state ODE*

- `run_gust_limited.jl` — 17 m/s Hann gust with limiter; β transient response
- `run_turbulent_limited.jl` — 15% TI at rated wind with limiter active

Demonstrates limiter robustness under realistic stochastic and transient loads.

### Priority 3 — Power limiter validation
*Engineering credibility*

- Plot β(t) vs P(t) time series from `run_wind_ramp_limited.jl`
- Verify settling time and steady-state β match design doc §2 table
- Check kp_elev tuning: currently 5×10⁻⁵ gives ~0.5°/s per kW overpower; β settles within 60s at 23 m/s (confirmed by `β range: 23.0°–52.5°`, `Mean last 30 s: 10.01 kW`)

### Priority 4 — 50 kW limiter validation
*Scale-up check*

- Run `run_wind_ramp_limited.jl` scaled to 50 kW (mass_scale applied)
- Verify kp_elev scaling (1/power_ratio) gives same relative response
- Check β_max = 67° is not exceeded

### Priority 5 — Force colouring contrast fix
*Visual bug — pre-tension collapses colour range*

Colour segments relative to their SWL (3500 N for tethers) rather than T_max_run. This ensures the full blue→red range is always visible regardless of the operating condition.

### Priority 6 — Extended telemetry CSV
*Enables post-run analysis*

Add per-segment tension columns T_seg_1..N and compression C_ring_1..M to CSV. force_analysis.jl already computes these — just needs CSV writer update in main.jl.

---

## File Map

```
TRPTKiteTurbineJulia/
├── src/
│   ├── parameters.jl       — SystemParams (24 fields), params_10kw/50kw, mass_scale
│   ├── aerodynamics.jl     — BEM Cp(λ)/CT(λ) tables, cp_at_tsr, ct_at_tsr
│   ├── wind_profile.jl     — Hellmann law, hub_altitude, wind constructors
│   ├── dynamics.jl         — 2-state ODEs, 3-state limiter ODEs, instantaneous_power
│   ├── geometry.jl         — compute_trpt_geometry, compute_blade_geometry, world_ground_plane
│   ├── force_analysis.jl   — ForceState, element_forces, run_force_scan
│   └── visualization.jl    — build_trpt_scene, _build_3d_axes!, _build_hud!, _build_controls!
├── scripts/                — 9 runnable simulation/recording/power-curve scripts
├── test/                   — 395 unit tests (all passing)
├── docs/plans/             — design docs, implementation plans, state snapshots
└── *.pdf / *.xlsx          — Source physics references (DRR, AeroDyn BEM, Mass Scaling)
```
