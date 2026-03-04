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
# → 168/168 tests pass
```

---

## Scripts

All scripts are run from the project root with `julia --project=. scripts/<name>.jl`.

### Interactive dashboards (GLMakie 3D window)

| Script | Description |
|--------|-------------|
| `run_dashboard_10kw.jl` | 10 kW prototype at rated wind (11 m/s) |
| `run_dashboard_50kw.jl` | 50 kW scaled configuration |
| `run_wind_ramp.jl` | Wind ramps from 6 → 14 m/s over t = 20–100 s |
| `run_gust.jl` | Smooth Hann-window gust peaking at 17 m/s at t = 60 s |
| `run_turbulent.jl` | IEC 61400-1 Class A turbulence, 15 % TI |

Press **▶ Play** in the sidebar, or drag the time slider manually. Press **Enter** in the terminal to exit.

### Video recording (headless, saves MP4)

| Script | Output |
|--------|--------|
| `record_simulation.jl` | `output/simulation_dashboard.mp4` — 10 kW, 120 s |
| `record_50kw.jl` | `output/simulation_50kw.mp4` — 50 kW, 120 s |

Uses the bundled `FFMPEG_jll` — no system ffmpeg needed.

### Power curve (headless, prints table + saves CSV)

```bash
julia --project=. scripts/power_curve.jl
# → output/power_curve.csv
```

Sweeps 4–16 m/s and reports steady-state ground power for both 10 kW and 50 kW configs.

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
  parameters.jl     — SystemParams struct, params_10kw(), params_50kw(), mass_scale()
  wind_profile.jl   — Hellmann wind shear + wind model constructors
                      (steady_wind, wind_ramp, gust_event, turbulent_wind)
  dynamics.jl       — trpt_ode!  (steady) and trpt_ode_wind! (time-varying wind)
  visualization.jl  — GLMakie 3D scene + telemetry HUD
  main.jl           — orchestration entry point

scripts/            — runnable simulation scenarios (see above)
test/               — 168 unit tests covering all source modules
docs/plans/         — original design documents
output/             — generated artefacts (gitignored except CSV)
```

---

## Key physics

- **Rotor**: ring of blades rotating at the top of the TRPT shaft
- **TRPT shaft**: tapered tensile tether; torsional stiffness k_i ∝ r_i (series springs)
- **Ground station**: PTO damper τ = c_pto × ω_ground
- **Wind shear**: Hellmann power law, exponent α = 1/7
- **Aerodynamic power**: P = ½ ρ v³ π R² Cp cos³(β), Cp = 0.22 (AeroDyn BEM)
- **Mass scaling**: P ∝ R² → R ∝ P^(1/2); c_pto ∝ P²

10 kW prototype at 11 m/s rated wind: **~9.1 kW mean**, **~9.6 kW peak**.
50 kW scaled config: **~45.6 kW mean**, **~48.8 kW peak**.
