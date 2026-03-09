# TRPT Kite Turbine Simulator

Julia simulation of a Tensile Rotary Power Transmission (TRPT) airborne wind energy system — a ring-rotor kite turbine that transfers torque to the ground via a tensile tether shaft.

Built by [Windswept & Interesting Ltd](https://windswept.energy) as part of the 10→50 kW kite turbine development programme.

---

## Requirements

- Julia 1.12+
- A display server for interactive GLMakie windows (X11 / Wayland)

## Installation

```bash
git clone https://github.com/rodWindswept/TRPTKiteTurbineJulia.git
cd TRPTKiteTurbineJulia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Run the tests

```bash
julia --project=. test/runtests.jl
# → 395/395 tests pass
```

---

## Scripts

All scripts are run from the project root with `julia --project=. scripts/<name>.jl`.

### Interactive dashboards (GLMakie 3D window)

| Script | Description |
|--------|-------------|
| `run_dashboard_10kw.jl` | 10 kW prototype at rated wind (11 m/s), fixed elevation β = 30° |
| `run_dashboard_50kw.jl` | 50 kW scaled configuration |
| `run_wind_ramp.jl` | Wind ramps 6 → 23 m/s over t = 20–100 s — 2-state MPPT, no limiter |
| `run_wind_ramp_limited.jl` | Wind ramps 6 → 23 m/s — **elevation power limiter active**: β rises from 23° as power exceeds rated |
| `run_gust.jl` | Smooth Hann-window gust peaking at 17 m/s at t = 60 s |
| `run_turbulent.jl` | IEC 61400-1 Class A turbulence, 15 % TI at 11 m/s mean |

Press **▶ Play** in the sidebar, or drag the time slider manually. Press **Enter** in the terminal to exit.

The **Live Telemetry** panel shows real-time elevation β — watch it rise above 23° in `run_wind_ramp_limited.jl` as the limiter engages above rated wind. The 3D shaft tilts accordingly.

The dashboard **Re-run** panel (unlock the toggle) lets you change wind speed, elevation β, and MPPT gain k_mppt on-the-fly and re-solve the 120 s ODE without restarting the script.

### Video recording (headless, saves MP4)

| Script | Output |
|--------|--------|
| `record_simulation.jl` | `output/simulation_dashboard.mp4` — 10 kW, 120 s |
| `record_50kw.jl` | `output/simulation_50kw.mp4` — 50 kW, 120 s |

Uses the bundled `FFMPEG_jll` — no system ffmpeg needed.

### Power curve (headless, prints table + saves CSV + PNG)

```bash
julia --project=. scripts/power_curve.jl
# → output/power_curve.csv  (4 curves: 10 kW MPPT, 10 kW limited, 50 kW MPPT, 50 kW limited)
# → output/power_curve.png
```

Sweeps 4–23 m/s and reports steady-state ground power for both configurations, in both MPPT-only and elevation-limiter modes. The limiter curves show a flat cap at rated power above rated wind.

### Full orchestration

```bash
julia --project=. src/main.jl        # 10 kW
julia --project=. src/main.jl 50     # 50 kW (or any rated kW)
```

Runs ODE, exports `output/telemetry_summary.csv`, prints a flight summary, and opens the interactive dashboard.

---

## Project structure

```
src/
  parameters.jl     — SystemParams struct (24 fields), params_10kw(), params_50kw(), mass_scale()
                      Includes 5 elevation-limiter fields: p_rated_w, β_min, β_max, β_rate_max, kp_elev
  aerodynamics.jl   — BEM Cp(λ)/CT(λ) lookup tables; cp_at_tsr(λ), ct_at_tsr(λ) interpolation
  wind_profile.jl   — Hellmann wind shear + wind model constructors
                      (steady_wind, wind_ramp, gust_event, turbulent_wind)
  dynamics.jl       — 2-state ODEs: trpt_ode!, trpt_ode_wind!, instantaneous_power
                      3-state ODEs: trpt_ode_limited!, trpt_ode_wind_limited! (elevation power limiter)
                      instantaneous_power/3 (3-arg variant for limiter trajectories)
  geometry.jl       — compute_trpt_geometry, compute_blade_geometry, world_ground_plane
  force_analysis.jl — ForceState, element_forces, run_force_scan
  visualization.jl  — GLMakie 3D scene + telemetry HUD + dynamic re-run panel
  main.jl           — orchestration entry point

scripts/            — 9 runnable simulation/recording/power-curve scripts
test/               — 395 unit tests covering all source modules
docs/plans/         — design documents and implementation plans
output/             — generated artefacts (gitignored except CSV)
```

---

## Key physics

- **Rotor**: ring of blades rotating at the top of the TRPT shaft
- **TRPT shaft**: tapered tensile tether; torsional stiffness k_i ∝ r_i (series springs)
- **Ground station**: MPPT quadratic torque law τ = k_mppt × ω² (eliminates bistability vs old linear c_pto)
- **Wind shear**: Hellmann power law, exponent α = 1/7
- **Aerodynamics**: TSR-dependent Cp(λ) and CT(λ) from AeroDyn BEM tables (NACA4412, 3-blade)
- **Aerodynamic power**: P = ½ ρ v³ π R² Cp(λ) cos³(β), Cp peak ≈ 0.22 at λ_opt ≈ 4.1
- **Elevation power limiter**: proportional controller raises β when P > P_rated; dβ/dt = clamp(kp_elev × P_error, bounded by β_rate_max); β clamped to [β_min=23°, β_max=67°]
- **Mass scaling**: P ∝ R² → R ∝ P^(1/2); k_mppt ∝ P^(5/2)

### 10 kW prototype (11 m/s rated wind, β = 30°, k_mppt = 11 N·m·s²/rad²)
- Steady-state: **~9.1 kW mean**, **~9.6 kW peak**
- Elevation limiter: β rises from 23° to ~52° at 23 m/s; mean power last 30 s ≈ 10.0 kW

### 50 kW scaled configuration
- Rotor radius 11.18 m, hub altitude 33.5 m
- Steady-state: **~45.6 kW mean**, **~48.8 kW peak**

---

## ODE state vectors

| Mode | State u | Description |
|------|---------|-------------|
| 2-state (MPPT) | `[α_tot, ω]` | TRPT twist angle + airborne rotor speed |
| 3-state (limiter) | `[α_tot, ω, β]` | As above + dynamic elevation angle |

The 3-state ODEs (`trpt_ode_limited!`, `trpt_ode_wind_limited!`) replace `p.elevation_angle` with the live state `u[3]` everywhere in the physics. All 2-state functions are unchanged.
