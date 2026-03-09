# Elevation-Based Power Limiter — Design Document

**Date:** 2026-03-09
**Status:** Approved — ready for implementation
**Author:** Rod Read / Claude Code

---

## 1. Motivation

The TRPT Kite Turbine simulator currently has no upper power bound. With MPPT tracking, power
grows as v³ above rated wind — the 10 kW config already outputs >20 kW at 15 m/s. In reality,
uncontrolled runaway power risks structural failure of the tether, ground station, and blades.

The chosen control strategy mirrors what would be implemented on real hardware: **raise the
elevation angle β** to increase cosine losses in the rotor disc. Since P_aero ∝ cos³(β), a
modest rise in β produces strong power shedding. This is a physically realisable action —
either by the lifter kite reducing lift (descending the rotor) or by relaxing mooring tension.

---

## 2. Physics Basis

Aerodynamic power in the current model (dynamics.jl step 6):

```
P_aero = 0.5 × ρ × v_hub³ × π × R² × Cp(λ) × cos³(β)
```

`cos³(β)` is already in the model. Raising β directly reduces aerodynamic torque without
changing TSR or Cp — the MPPT continues to track λ_opt, so rotor efficiency is preserved
while total power intake is reduced.

Power shedding headroom (10 kW config, v_rated = 11 m/s):

| β     | cos³(β) | P_aero relative |
|-------|---------|-----------------|
| 23°   | 0.778   | 100% (minimum, rated baseline) |
| 35°   | 0.551   | 71% |
| 45°   | 0.354   | 46% |
| 55°   | 0.188   | 24% |
| 67°   | 0.059   | 8%  |

At v = 20 m/s (the top of the ramp scenario), P_aero_unlimited ≈ (20/11)³ × P_rated ≈ 6×
rated. β ≈ 55° would restore rated power at this extreme.

---

## 3. State Vector Change

The ODE grows from 2 to 3 states:

| Index | Variable | Unit | Description |
|-------|----------|------|-------------|
| `u[1]` | α_tot | rad | TRPT total twist angle — unchanged |
| `u[2]` | ω | rad/s | Airborne rotor angular velocity — unchanged |
| `u[3]` | β | rad | Elevation angle — **new dynamic state** |

Everywhere `p.elevation_angle` appears inside the ODE physics
(`hub_altitude`, `cos(β)^3` in P_aero, `sin(β)` in tether drag) it is replaced by `u[3]`.

Initial condition: `[0.0, ω₀, β_min]` — starts at minimum safe elevation.

`p.elevation_angle` is retained in `SystemParams` as the **nominal/reference** elevation used
by:
- The existing 2-state ODEs (unchanged)
- The dashboard β slider initial value in **manual re-run mode**
- When the dashboard is in **auto limiter mode**, β₀ = p.β_min and evolves dynamically

---

## 4. SystemParams — New Fields

```julia
p_rated_w  ::Float64   # Rated electrical power (W). Default: 10_000.0 (10 kW)
β_min      ::Float64   # Minimum elevation angle (rad). Default: deg2rad(23.0)
                       #   Safety floor: blade tips must clear the ground
β_max      ::Float64   # Maximum elevation angle (rad). Default: deg2rad(67.0)
                       #   Lifter kite pull limit (IRL varies with wind speed)
β_rate_max ::Float64   # Maximum elevation change rate (rad/s). Default: deg2rad(1.0)
                       #   1°/s — representative of lifter kite / mooring servo speed
kp_elev    ::Float64   # Proportional gain (rad/W/s). Default: 5e-5
                       #   Gives ≈0.5°/s response per 1 kW overpower
```

`mass_scale` scales `p_rated_w` linearly with power ratio. All angle fields are fixed
(independent of scale). `kp_elev` scales as `power_ratio^(-1)` so fractional overpower gives
the same angular response rate regardless of turbine size.

---

## 5. Control Law

```
P_gen   = τ_transmitted × ω_ground          (computed inline in ODE — already available)
P_error = P_gen − p.p_rated_w               (positive → over-rated → raise β)

rate_raw = p.kp_elev × P_error

lower_rate = (u[3] ≤ p.β_min) ? 0.0 : −p.β_rate_max
upper_rate = (u[3] ≥ p.β_max) ? 0.0 : +p.β_rate_max

du[3] = clamp(rate_raw, lower_rate, upper_rate)
```

Below rated power P_error < 0 → β decreases back toward β_min (MPPT region recovery).
Above rated power P_error > 0 → β increases to shed load.
Saturation guards prevent β winding outside [β_min, β_max] regardless of solver step size.

No integrator term in the first implementation. Steady-state β will settle slightly above the
ideal β_rated — conservative (slightly under-rated), safe, and simple to reason about.

---

## 6. New ODE Functions

Two new functions added to `src/dynamics.jl`, keeping all existing 2-state functions intact:

```julia
trpt_ode_limited!(du, u, p, t)           # 3-state, steady wind
trpt_ode_wind_limited!(du, u, pw, t)     # 3-state, time-varying wind
```

The existing `trpt_ode!`, `trpt_ode_wind!`, and `instantaneous_power` are **not modified**.
All existing tests, scripts, and dashboards continue to work without change.

`instantaneous_power` is extended with an optional `β` argument:

```julia
instantaneous_power(p, u)          # existing 2-state: uses p.elevation_angle
instantaneous_power(p, u, β)       # new 3-state variant: uses β = u[3]
```

---

## 7. Dashboard Integration

### Manual mode (existing behaviour)
- β slider sets the initial `p.elevation_angle` used in the 2-state ODE re-run
- β slider range corrected: minimum 23° (was incorrectly allowing lower values)
- No change to existing HUD layout

### Auto limiter mode (new)
- A **"Auto limiter" toggle** is added to the re-run panel
- When enabled, re-run uses `trpt_ode_limited!` (3-state)
- Initial β = β_min (23°), not the slider value
- β slider becomes read-only display of current β (greyed out)
- HUD gains a fourth panel row: `β (°)` showing the elevation time series

The β(t) trace in the 3D scene: the inclined shaft angle updates with `u[3]` during playback,
so the viewer watches the kite rise and fall as the limiter acts.

---

## 8. New Script

`scripts/run_wind_ramp_limited.jl`

- Wind: 6 → 20 m/s linear ramp over 120 s
- ODE: `trpt_ode_wind_limited!`, 3-state
- Warm-start: `ω₀ = 4.1 × v_hub(6 m/s) / R`, `β₀ = deg2rad(23.0)`
- Output: interactive GLMakie dashboard via `build_trpt_scene`
- Trajectory carries `beta = sol[3, :]` alongside `alpha_tot`, `omega`, `power_kw`, `v_hub`
- The β trace is plotted in the HUD torque panel (or a dedicated row) as "Elevation β (°)"

Expected behaviour: flat 10 kW output above ~11 m/s; β rises from 23° toward ~55° at 20 m/s.

---

## 9. Testing

New tests in `test/test_dynamics.jl`:

- **Limiter active at high wind**: solve 60 s at v = 20 m/s, assert `mean(P_gen[end-30s:end]) ≤ p.p_rated_w × 1.05`
- **β stays in bounds**: assert `all(β_trace .≥ p.β_min - 1e-6)` and `all(β_trace .≤ p.β_max + 1e-6)`
- **MPPT region unaffected**: at v = 8 m/s (below rated), assert β remains at β_min throughout
- **Recovery after gust**: after a 20 m/s gust, β returns toward β_min once wind drops below rated

---

## 10. Out of Scope (Future)

- **Integral term** — steady-state β error is small and conservative; PI controller is a
  natural next step once proportional behaviour is validated
- **Blade pitch / Cp reduction** — physically richer but requires BEM table extension
- **Lifter wind-speed dependency on β_max** — IRL β_max varies with wind speed as lifter
  lift changes; modelling this requires a lifter kite aerodynamic sub-model
- **Cut-out control** — above some extreme wind speed, full shutdown (β → 90°, TRPT released)
