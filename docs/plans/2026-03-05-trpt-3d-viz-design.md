# TRPT Kite Turbine — High-Fidelity 3D Visualization Engine Design

**Date:** 2026-03-05
**Author:** Rod Read / Claude Code
**Status:** Approved

---

## 1. Purpose

Extend the existing Julia/GLMakie simulation with a high-fidelity, interactive 3D validation
environment for engineers to verify kinematic integrity, spatial alignment, and structural load
distribution of the TRPT Kite Turbine system across its full range of elevation angles and
dynamic torque conditions.

---

## 2. Architecture — Layered Module Refactor (Option C)

```
src/
├── parameters.jl          # unchanged
├── wind_profile.jl        # unchanged
├── dynamics.jl            # unchanged (ODE, instantaneous_power)
├── geometry.jl            # NEW — 3D geometry, replaces geometry portion of visualization.jl
│                          #   compute_trpt_geometry(p, alpha_tot, shaft_dir)
│                          #   compute_blade_geometry(p, alpha_tot, shaft_dir)
│                          #   world_ground_plane()
├── force_analysis.jl      # NEW — per-element forces from ODE state
│                          #   element_forces(p, u, v_hub) -> ForceState
├── visualization.jl       # REFACTORED — rendering & UI only
│                          #   build_trpt_scene(p, traj) — enhanced
└── main.jl                # unchanged (orchestration)

test/
├── test_geometry.jl       # expanded — inclined shaft, blades
├── test_force_analysis.jl # NEW
└── ... (others unchanged)
```

### Coordinate system (shared across all modules)

- `+x` = downwind (wind direction projected onto ground plane)
- `+y` = crosswind right
- `+z` = vertical up
- TRPT ground anchor at world origin `(0, 0, 0)`
- `shaft_dir = [cos(φ_w)·cos(β), sin(φ_w)·cos(β), sin(β)]` — unit vector along shaft axis
  where `φ_w` = wind azimuth, `β` = elevation angle

---

## 3. Geometry Module (`src/geometry.jl`)

### `compute_trpt_geometry(p, alpha_tot, shaft_dir)`

**Current problem:** every ring is centred at `(0, 0, z)` — the shaft appears vertical
regardless of elevation angle.

**Fix:** each ring centre steps along `shaft_dir` by `l_seg` per level. Ring node positions
are offset radially in the plane perpendicular to `shaft_dir`:

```
ring_centre[i] = (i-1) × l_seg × shaft_dir
perp_axes      = two orthonormal vectors ⊥ shaft_dir (via cross products, computed once)
node[i,j]      = ring_centre[i] + r_i × (cos(φ_ij) × perp_axes[1]
                                        + sin(φ_ij) × perp_axes[2])
```

where `φ_ij = theta_i + (j-1) × 2π/n_lines`, taper and twist logic unchanged.
Return shape stays `(n_levels, n_lines, 3)` — all existing geometry tests hold.

### `compute_blade_geometry(p, alpha_tot, shaft_dir)`

Returns `(n_blades, 4, 3)` — four corners of each rectangular blade quad.

**Blade extent** (key correction — ring is NOT at blade root):

```julia
const BLADE_INNER_FRAC = 0.30   # ~30% of blade span is inboard of the ring

blade_span     = (p.rotor_radius - p.trpt_hub_radius) / (1.0 - BLADE_INNER_FRAC)
blade_inner_r  = p.trpt_hub_radius - BLADE_INNER_FRAC * blade_span   # inboard tip
blade_outer_r  = p.rotor_radius                                        # outboard tip
chord          = p.rotor_radius * 0.15                                 # 15% radius
```

For 10 kW: `blade_span ≈ 4.29 m`, `blade_inner_r ≈ 0.71 m`, ring at `trpt_hub_radius = 2.0 m`.

Blade count: `p.n_blades` (design principle: `n_blades = n_lines = polygon sides` for even
loading — the 10 kW prototype has `n_blades = 3`, `n_lines = 5` as separate tuneable params).

Each blade is locked in the plane perpendicular to `shaft_dir` (rotor plane), rotated by
`alpha_tot + (b-1) × 2π/n_blades` — satisfies the kinematic constraint at all elevation angles.

### `world_ground_plane()`

Grid of lines at `z = 0`, ±40 m in x and y. Rendered grey, semi-transparent.

---

## 4. Force Analysis Module (`src/force_analysis.jl`)

### `ForceState` struct

```julia
struct ForceState
    tether_tension   :: Vector{Float64}   # n_seg values — axial tension per segment (N)
    ring_compression :: Vector{Float64}   # n_rings values — hoop compression per ring (N)
    tau_transmitted  :: Float64           # total TRPT torque (N·m)
    tau_aero         :: Float64           # aerodynamic torque (N·m)
    tau_drag         :: Float64           # tether drag torque (N·m)
end
```

### `element_forces(p, u, v_hub) -> ForceState`

Derived analytically from ODE state `u = [alpha_tot, omega]`:

**Tether tension per segment:**
```
T_i = (k_i × sin(α_i)) / r_i   +   m_seg × ω² × r_i   +   m_above_i × g × sin(β)
```
- Torsional component: torque arm at radius `r_i`
- Centrifugal: tether segment mass × ω² × r_i
- Gravitational: mass above segment projected along shaft axis

**Ring hoop compression** (polygon geometry):
```
C_i = T_i × r_i / (2 × sin(π / n_lines))
```

### Colour scale — auto-calibrated per run with safety margin

Since the full ODE trajectory is pre-computed before playback, scan all frames after each
solve to find `T_max_run` and `C_max_run`. Scale set once per solve.

**Colourbar layout** (fixed in scene, always visible):

```
 SWL ────────┤                        ├── SWL
             │  ░░░░░░░░░░░░░░░░░░░░  │   grey extension: max_run → SWL
             │  FoS = SWL / max_run   │   factor of safety label
 max_run ────┤                        ├── max_run
             │  ████████████████████  │   blue → red: 0 to max_run
 0 ──────────┤                        ├── 0
```

- Blue→red gradient spans `0 → max_run` — full dynamic range of observed loads
- Grey extension shows remaining headroom to SWL
- FoS number displayed on the bar
- If `max_run > SWL`: bar inverts, label turns red (overload warning)

**Reference SWL values:**
- Tether: Dyneema 3 mm ≈ 3500 N (DRR §5.2)
- Ring strut: `RING_SWL = 500.0` N (CFRP tube buckling, conservative — adjustable constant)

---

## 5. Visualization & UI (`src/visualization.jl`)

### Layout

```
┌─────────────────────────────────┬──────────────────────────────────┐
│                                 │  TELEMETRY HUD                   │
│       3D VIEWPORT               │  t, P, ω, α_seg, collapse margin │
│       (Axis3, interactive)      │  V_hub, wind azimuth             │
│                                 ├──────────────────────────────────┤
│  - inclined TRPT shaft          │  FORCE COLOURBAR                 │
│  - tether lines (force-coloured)│  tether tension + ring compr.    │
│  - ring polygons (force-coloured│  with FoS extension              │
│  - n_blades rectangular blades  ├──────────────────────────────────┤
│  - ground plane grid            │  CONTROLS                        │
│  - wind direction arrow         │  [ elevation β   0 ──── 75° ]   │
│                                 │  [ wind azimuth  0 ─── 360° ]   │
│                                 │  [ time slider              ]   │
│                                 │  [ ▶ Play ]  [ solver ▾ ]       │
│                                 ├──────────────────────────────────┤
│                                 │  DYNAMIC TORQUE MODE             │
│                                 │  [ ☐ Enable ]                   │
│                                 │  [ c_pto  (N·m·s/rad) slider ]  │
│                                 │  [ Re-run ODE ]                  │
└─────────────────────────────────┴──────────────────────────────────┘
```

### Reactive behaviours

**Elevation & azimuth sliders** — recompute `shaft_dir`, update geometry observables
instantly. No ODE re-run (purely kinematic transform).

**Time slider / Play** — drives frame index observable, unchanged from current behaviour.

**Solver dropdown** (`Menu`): `Tsit5` (default), `RK4`, `Euler`. Stored as `Observable{Symbol}`.

**Dynamic Torque panel:**
1. Checkbox enables the panel
2. `c_pto` slider: 100 → 50 000 N·m·s/rad (log scale), default = `p.c_pto`
3. "Re-run ODE": rebuild `SystemParams` with new `c_pto`, solve with chosen solver,
   rescan forces → recalibrate colourbar, update all observables atomically, reset to frame 1

### Internal structure

`build_trpt_scene()` delegates to three helpers:
- `_build_3d_axes!(fig, p, ...)` — Axis3, geometry renderables, ground plane, wind arrow
- `_build_hud!(fig, ...)` — telemetry labels, colourbar
- `_build_controls!(fig, ...)` — sliders, buttons, dynamic torque panel

---

## 6. Data Synchronisation Schema

Extended telemetry CSV (appended columns, backward-compatible):

```
t, alpha_tot, omega, power_kw, v_hub,
T_seg_1 .. T_seg_N,          # tether tension per segment (N), N = n_rings+1
C_ring_1 .. C_ring_M,        # ring hoop compression (N), M = n_rings
tau_aero, tau_drag, tau_transmitted,
T_max_run, C_max_run         # run-calibrated scale maxima (constant per run)
```

**Verifiability:** at any frame `i`, 3D node positions are fully determined by `alpha_tot[i]`,
`elevation_angle`, and `wind_azimuth`. Element colours are fully determined by
`T_seg_j[i] / T_max_run` and `C_ring_k[i] / C_max_run`. An engineer can reproduce any frame
from the CSV alone.

---

## 7. Future Work (out of scope for this build)

- **Parametric geometry screen** — pre-run UI panel to adjust `n_lines`, `n_rings`,
  `tether_length`, `rotor_radius`, etc. before solving
- **Variable elevation angle** — elevation responding dynamically to wind load / kite lift
- **Terrain elevation** — offset TRPT ground point from a DEM for site siting visualisations
