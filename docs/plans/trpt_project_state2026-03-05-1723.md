# TRPTKiteTurbineJulia — State of Project & Development Priorities

*Analysis date: 2026-03-05*

---

## What Has Been Built

### Phase 1 — Core Simulator ✅ Complete
Documented in [2026-03-03-trpt-simulator-design.md](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/docs/plans/2026-03-03-trpt-simulator-design.md) and its implementation plan — all tasks done.

| Module | File | Status |
|--------|------|--------|
| Physical parameters + scaling | [src/parameters.jl](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/src/parameters.jl) | ✅ Implemented & evolved |
| Hellmann wind profile | [src/wind_profile.jl](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/src/wind_profile.jl) | ✅ Complete |
| TRPT ODE dynamics | [src/dynamics.jl](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/src/dynamics.jl) | ✅ Complete |
| Entry point / orchestration | [src/main.jl](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/src/main.jl) | ✅ Complete |
| Test suite | `test/` (168 tests) | ✅ All pass |
| Simulation scripts | `scripts/` (10 scripts) | ✅ All runnable |

### Phase 2 — High-Fidelity 3D Viz Engine ✅ Complete
Documented in [2026-03-05-trpt-3d-viz-design.md](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/docs/plans/2026-03-05-trpt-3d-viz-design.md) and implementation plan — all 3 tasks done.

| Module | File | Status |
|--------|------|--------|
| Inclined shaft 3D geometry + blades | [src/geometry.jl](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/src/geometry.jl) | ✅ Implemented |
| Per-element force analysis (ForceState) | [src/force_analysis.jl](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/src/force_analysis.jl) | ✅ Implemented |
| Refactored GLMakie scene + UI | [src/visualization.jl](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/src/visualization.jl) | ✅ Implemented |

**Delivered features:**
- Inclined shaft geometry (ring planes ⊥ shaft axis, correct elevation)
- Tapered TRPT (r_bottom → r_top taper derived from `trpt_rL_ratio`)
- Rectangular rotor blades (5 blades, `BLADE_INNER_FRAC = 0.30` geometry)
- Force-coloured tether lines (blue→red) and ring polygons — **code present but colour gradient not visually working** (likely lifter pre-tension dominates all segments, collapsing colour contrast)
- FoS HUD text labels — **implemented but not confirming visible in dashboard** (may need a frame tick to trigger; possible layout scroll issue)
- Elevation β slider (0–75°), wind azimuth φ slider (0–360°) — kinematic only
- Dynamic torque panel: c_pto log-scale slider + ODE re-run on demand
- Solver dropdown: Tsit5 / RK4 / Euler — **no in-UI explanation of solver differences** (UX gap)
- Lifter kite force arrows at top ring nodes (gold lines)
- Play/Pause + time scrub

---

## Key Parameter Changes Since Original Design

These are important — the implementation diverged from the original design document in physics-informed ways:

| Parameter | Original Design | Current Implementation | Reason |
|-----------|-----------------|----------------------|--------|
| `tether_length` | 150 m | **30 m** | DRR §5.2 — actual hardware is 30 m TRPT |
| `n_lines` | 6 | **5** | DRR §5.2 — 5 Dyneema tethers |
| `tether_diameter` | 4 mm | **3 mm** | DRR §5.2 — Dyneema type 01505 |
| `n_rings` | 10 | **14** | 30 m / 2 m spacing − 1 = 14 (DRR) |
| `m_ring` | 2.5 kg/ring | **0.4 kg/ring (flat)** | DRR §5.2 hardware measurement — but should be physics-derived from ring radius and hoop compression capacity (see Priority 6 below) |
| `n_blades` | 3 | **5** | One blade per tether line (n_blades = n_lines) |
| `cp` | 0.15 | **0.22** | AeroDyn BEM (NACA4412 3-blade); replaces proxy |
| `h_ref` | 150 m | **15 m** | Hub altitude = 30 × sin(30°) = 15 m |
| Geom. scaling exponent | (P/P₀)^(1/3) | **(P/P₀)^(1/2)** | P ∝ R² (swept area), not R³ |
| New: `trpt_hub_radius` | — | **2.0 m** | Rotor-end ring radius (DRR Grasshopper FEA) |
| New: `trpt_rL_ratio` | — | **0.74** | r/L geometry constraint per segment |
| New: `lifter_elevation` | — | **70°** | Lifter kite line angle |

> [!IMPORTANT]
> The `h_ref = 15 m` at 11 m/s rated wind (not 10 m/s at 150 m) is a significant change. The hub is low — Hellmann shear over 15 m has very little effect. The 150 m original was likely a modelling fiction.

---

## Outstanding Work — The Todo Note

From [docs/plans/todo 1 implement colour.txt](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/docs/plans/todo%201%20implement%20colour.txt):

1. **Rotor power dynamics validation** — verify the power curve is realistic vs wind speed and elevation angle (per prior research / AeroDyn BEM data)
2. **Colour change for compression range guide** — implement the FoS colourbar visual as specified in the viz design (grey extension above `max_run` up to SWL, with FoS number)
3. **Rotary solution dynamics** — show rotor spin speed and TRPT twist rate for each combination of generator torque (`c_pto`) and wind input

---

## Gaps vs Approved Design Docs

### From [2026-03-05-trpt-3d-viz-design.md](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/docs/plans/2026-03-05-trpt-3d-viz-design.md) — Not Yet Implemented

| Feature | Design Section | Status |
|--------|---------------|--------|
| Force colour gradient visible in viewport | §4 | ❌ Code present; not working visually — pre-tension collapses contrast |
| FoS HUD text labels updating | §5 | ⚠️ Code present; not confirmed working in dashboard |
| In-UI solver explanation (Tsit5/RK4/Euler) | §5 UI | ❌ No explanatory text next to dropdown |
| Visual colourbar widget (gradient bar in scene) | §4 & §5 | ❌ Text-only — no graphical bar |
| Grey extension above `max_run` → SWL on bar | §4 colour scale | ❌ Not done (todo item 2) |
| FoS number on graphical colourbar | §4 colour scale | ❌ Not done |
| Extended telemetry CSV (`T_seg_i`, `C_ring_j`, ...) | §6 | ❌ Current CSV is basic |
| Parametric geometry screen (pre-run UI) | §7 Future work | Explicitly deferred |
| Variable elevation responding to wind load | §7 Future work | Explicitly deferred |

### Physics Model Known Limitations (documented)

From design doc §6:
1. **Rigid body** — no blade aeroelastic deformation
2. **Bulk Cp = 0.22** — full BEM loop not implemented per-frame
3. **Lift kite decoupled** — TRPT axis assumed static (lifter is visualised but not dynamically coupled)

---

## Suggested Development Priorities for Discussion

### Priority 1 — Power Dynamics Validation (Todo #1)
*Low risk, high value for commercial credibility*

- Generate and plot the power curve from [scripts/power_curve.jl](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/scripts/power_curve.jl) output
- Cross-check P vs V_wind against AeroDyn BEM data in [Rotor_TRTP_Sizing_Iteration2.xlsx](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/Rotor_TRTP_Sizing_Iteration2.xlsx)
- Verify the cos³(β) elevation correction matches expectations at 20°, 30°, 45°
- Check the rated-wind operating point: at 11 m/s, 30°, c_pto=5000 → expected ≈ 9–10 kW

### Priority 2 — Fix Force Colouring + HUD (Bugs)
*These are confirmed not working — fix before the graphical colourbar*

**Force colour contrast issue:** The lifter pre-tension `T_lifter = m_airborne×g/sin(70°)/n_lines` is applied equally to every segment, which dominates and eliminates per-segment contrast. Fix options:
  - Colour relative to SWL instead of `T_max_run` — always shows meaningful safety margin
  - Or subtract the constant pre-tension from the display range so only the variable component drives colour

**FoS HUD:** Confirm whether labels are visible; may need initial trigger on scene open

**Solver dropdown UX:** Add a static Label below the Menu explaining:
> Tsit5 – adaptive, accurate (default) · RK4 – fixed-step, predictable · Euler – basic, fast but may drift

### Priority 2b — Graphical Colourbar (Todo #2)
*Only meaningful once force colouring is working*

- Replace text force labels with visual gradient bar: blue→red (0→T_max_run), grey extension to SWL, FoS number
- Use GLMakie `Colorbar` or custom rectangles

### Priority 3 — Rotary Dynamics Display (Todo #3)
*Most analytically interesting for engineering insight*

- Steady-state ω vs (V_wind, c_pto) surface plot — "operating map"
- Show twist accumulation rate dα/dt as a function of loading
- Could be a separate script or a new panel in the dashboard

### Priority 4 — Extended Telemetry CSV
*Enables post-run analysis and reproducibility*

- Add per-segment tension columns `T_seg_1..N` and compression `C_ring_1..M` to CSV
- Backward-compatible extension; [force_analysis.jl](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/src/force_analysis.jl) already computes these

### Priority 5 — Physics Refinements (Longer-term)
- Lock `h_ref` to actual hub altitude consistently (currently done; verify all scripts)
- Consider whether wind shear model is worth improving (15 m hub → very flat profile)
- Tether drag: the Tulloch model is an approximation; validate against 3D solver if available
- Lift kite dynamic coupling — major physics extension, requires separate project scoping

### Priority 6 — Physics-Derived Ring Mass Model
*Blocked on aerodynamics being settled first — do after Priority 1*

Each ring's mass should come from the structural requirement, but the load direction is not simply hoop compression (see Physics Notes below). Defer until rotor aero and blade geometry are resolved.

---

## Physics Notes — Real Geometry Effects Not Yet Modelled

> [!IMPORTANT]
> These are important real-world effects flagged for future modelling. Not bugs — documented physics gaps.

### Blade coning / bank angle
Real blades are slightly **banked**: outer tips swept toward the generator/ground end, inner tips toward the sky. In rotation this produces a **coning geometry** — the swept disc is not flat but conical. Effects:
- Effective rotor radius under centripetal load is **larger** than the static `rotor_radius`
- The `cos³(β)` inflow correction in the current power formula is for shaft inclination only, not blade coning — the two interact
- AeroDyn BEM `Cp = 0.22` may embed a specific coning assumption; needs checking against the BEM input geometry

### Ring loading direction — not always hoop compression
[force_analysis.jl](file:///home/rod/Documents/GitHub/TRPTKiteTurbineJulia/src/force_analysis.jl) currently assumes tether lines pull outward, ring struts resist in **compression**. However:
- Centripetal acceleration at high ω pulls the top ring **outward** (wider), which can put ring struts in **tension** rather than compression — especially at high spin rates
- The blade bank angle reinforces this: centripetal force on banked blade mass has a component increasing the coning angle dynamically, further loading the ring outward
- If ring struts are in tension, material choice changes (CFRP tube is poor in tension vs. cable/rope)
- The `ring_compression` output in `ForceState` may be the wrong sign / wrong quantity at high ω

**Implication for ring mass model:** The CFRP-buckling derivation (Priority 6) assumes compression. If the ring is centrifugally loaded in tension at operating ω, the sizing calculation and material choice are different. Resolve aero + blade geometry first.

---

## File Map

```
TRPTKiteTurbineJulia/
├── src/
│   ├── parameters.jl       — SystemParams, params_10kw/50kw, mass_scale
│   ├── wind_profile.jl     — Hellmann law, hub_altitude, wind_at_hub
│   ├── dynamics.jl         — trpt_ode!, trpt_ode_wind!, instantaneous_power
│   ├── geometry.jl         — compute_trpt_geometry, compute_blade_geometry, world_ground_plane
│   ├── force_analysis.jl   — ForceState, element_forces, run_force_scan
│   └── visualization.jl    — build_trpt_scene, _build_3d_axes!, _build_hud!, _build_controls!
├── scripts/                — 10 runnable simulation/recording/power-curve scripts
├── test/                   — 168 unit tests (all passing)
├── docs/plans/             — 4 design docs + 1 todo note
└── *.pdf / *.xlsx          — Source physics references
```
