# Elevation Power Limiter Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a proportional elevation-angle power limiter as a 3rd ODE state so the simulator enforces rated power above rated wind without touching any existing 2-state ODE functions or tests.

**Architecture:** Two new 3-state ODE functions (`trpt_ode_limited!`, `trpt_ode_wind_limited!`) replace `p.elevation_angle` with `u[3]` and append `dβ/dt` to the derivative. Five new fields in `SystemParams` carry the limiter parameters. The dashboard gains an "Auto limiter" toggle that switches the re-run between 2-state (manual β) and 3-state (dynamic β) modes.

**Tech Stack:** Julia 1.12, DifferentialEquations.jl, GLMakie.jl, Test.jl

**Design doc:** `docs/plans/2026-03-09-elevation-power-limiter-design.md`

---

## Task 1: Add limiter fields to SystemParams

**Files:**
- Modify: `src/parameters.jl`

Fields to append to the end of the `SystemParams` struct (after `k_mppt`):

```julia
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
```

**Step 1: Add 5 fields to the struct**

Append after `k_mppt::Float64` in the struct definition:

```julia
    p_rated_w  ::Float64
    β_min      ::Float64
    β_max      ::Float64
    β_rate_max ::Float64
    kp_elev    ::Float64
```

**Step 2: Update `params_10kw()` — append 5 values**

After `11.0,  # k_mppt` in the `SystemParams(...)` call, append:

```julia
        10_000.0,        # p_rated_w (W) — 10 kW rated
        deg2rad(23.0),   # β_min — blade-tip ground clearance floor
        deg2rad(67.0),   # β_max — lifter kite ceiling
        deg2rad(1.0),    # β_rate_max (rad/s) — 1 °/s actuation rate limit
        5e-5,            # kp_elev (rad/W/s) — 0.5 °/s per 1 kW overpower
```

**Step 3: Update `mass_scale()` — append 5 scaling rules**

After `base.k_mppt * power_ratio^2.5,` in the `SystemParams(...)` call, append:

```julia
        base.p_rated_w  * power_ratio,        # rated power scales linearly
        base.β_min,                            # angle does not scale
        base.β_max,                            # angle does not scale
        base.β_rate_max,                       # rate does not scale
        base.kp_elev    / power_ratio,         # larger turbine: same °/s per fractional overpower
```

**Step 4: Run full test suite**

```bash
julia --project=. test/runtests.jl
```

Expected: **384 tests pass** (struct change is additive — positional constructors in test files use `params_10kw()` not raw struct literals).

**Step 5: Commit**

```bash
git add src/parameters.jl
git commit -m "feat: add 5 limiter fields to SystemParams (p_rated_w, β_min, β_max, β_rate_max, kp_elev)"
```

---

## Task 2: Fix elevation slider minimum to 23°

**Files:**
- Modify: `src/visualization.jl:292`

The current slider range `0.0:1.0:75.0` allows values below the safe floor. Fix to start at 23°.

**Step 1: Fix the slider range**

```julia
# Before (line 292):
elev_slider = Slider(layout[3, 1]; range=0.0:1.0:75.0,
                     startvalue=rad2deg(p.elevation_angle))

# After:
elev_slider = Slider(layout[3, 1]; range=23.0:1.0:75.0,
                     startvalue=clamp(rad2deg(p.elevation_angle), 23.0, 75.0))
```

**Step 2: Run test suite**

```bash
julia --project=. test/runtests.jl
```

Expected: **384 tests pass**

**Step 3: Commit**

```bash
git add src/visualization.jl
git commit -m "fix: clamp elevation slider minimum to 23° (blade-tip clearance floor)"
```

---

## Task 3: Add 3-state limiter ODEs to dynamics.jl

**Files:**
- Modify: `src/dynamics.jl`

Add two new functions **after** `instantaneous_power`. Do NOT touch `trpt_ode!`, `trpt_ode_wind!`, or `instantaneous_power`.

The 3-state vector: `u[1]` = α_tot, `u[2]` = ω, `u[3]` = β (elevation, rad).
Everywhere `p.elevation_angle` is used inside the physics, use `u[3]` instead.

**Step 1: Write failing tests first**

Add to `test/test_dynamics.jl`:

```julia
@testset "trpt_ode_limited! — basic call (3-state)" begin
    p  = params_10kw()
    du = zeros(3)
    u  = [0.1, 1.0, deg2rad(23.0)]   # α, ω, β at minimum

    trpt_ode_limited!(du, u, p, 0.0)

    @test isfinite(du[1])
    @test isfinite(du[2])
    @test isfinite(du[3])
end

@testset "trpt_ode_limited! — β stays at floor under low wind" begin
    p  = params_10kw()
    # Low wind (4 m/s) → P_gen < P_rated → limiter tries to decrease β → clamps at β_min
    using DifferentialEquations
    sol = solve(ODEProblem(trpt_ode_limited!,
                           [0.0, 4.1*4.0/p.rotor_radius*0.8, p.β_min],
                           (0.0, 30.0), p),
                Tsit5(); reltol=1e-6, abstol=1e-6, saveat=1.0)
    β_trace = sol[3, :]
    @test all(β_trace .>= p.β_min - 1e-4)   # never goes below floor
end

@testset "trpt_ode_limited! — β rises under high wind" begin
    p  = params_10kw()
    # High wind (18 m/s): P_aero >> P_rated → β must rise above β_min
    using DifferentialEquations
    ω0 = 4.1 * 18.0 / p.rotor_radius
    sol = solve(ODEProblem(trpt_ode_limited!,
                           [0.0, ω0, p.β_min],
                           (0.0, 60.0), p),
                Tsit5(); reltol=1e-6, abstol=1e-6, saveat=1.0)
    β_trace = sol[3, :]
    @test maximum(β_trace) > p.β_min + deg2rad(2.0)   # β rises meaningfully
    @test all(β_trace .<= p.β_max + 1e-4)              # never exceeds ceiling
end

@testset "trpt_ode_limited! — power clamped near rated at high wind" begin
    p  = params_10kw()
    using DifferentialEquations, Statistics
    ω0 = 4.1 * 18.0 / p.rotor_radius
    sol = solve(ODEProblem(trpt_ode_limited!,
                           [0.0, ω0, p.β_min],
                           (0.0, 120.0), p),
                Tsit5(); reltol=1e-6, abstol=1e-6, saveat=1.0)
    # Average power over final 30 s must be within 15% of rated (limiter settling)
    idx = findall(t -> t >= 90.0, sol.t)
    P_vals = [instantaneous_power(p, [sol[1,i], sol[2,i]], sol[3,i]) for i in idx]
    @test mean(P_vals) <= p.p_rated_w * 1.15
end

@testset "instantaneous_power — 3-arg variant" begin
    p = params_10kw()
    # With explicit β, result must differ from default (p.elevation_angle = 30°)
    # when β is different
    P_default = instantaneous_power(p, [0.5, 2.0])
    P_high_β  = instantaneous_power(p, [0.5, 2.0], deg2rad(60.0))
    @test P_default isa Float64
    @test P_high_β  isa Float64
    @test P_high_β < P_default   # higher β → more cosine loss → less power
end
```

**Step 2: Run tests — expect failure**

```bash
julia --project=. test/runtests.jl
```

Expected: new tests FAIL with `UndefVarError: trpt_ode_limited!` and `MethodError: instantaneous_power`.

**Step 3: Extend `instantaneous_power` with optional β argument**

Locate `instantaneous_power` in `src/dynamics.jl`. Add a 3-argument method **directly after** the existing 2-argument function. The 3-argument method is identical but substitutes `β` for `p.elevation_angle` in the `P_aero` cos³ term:

```julia
"""
    instantaneous_power(p, u, β) -> Float64

3-state variant: as `instantaneous_power(p, u)` but uses elevation angle `β` (rad)
from the dynamic state `u[3]` instead of the fixed `p.elevation_angle`.
"""
function instantaneous_power(p::SystemParams, u::Vector{Float64}, β::Float64)::Float64
    n_segments   = p.n_rings + 1
    r_top        = p.trpt_hub_radius
    r_bottom     = 2.0 * p.tether_length * p.trpt_rL_ratio / n_segments - r_top
    spring_coeff = p.e_modulus * π * (p.tether_diameter / 2)^2 * p.n_lines * p.trpt_rL_ratio
    sum_inv_k    = sum(1.0 / (spring_coeff * (r_bottom + i / (n_segments - 1) * (r_top - r_bottom)))
                       for i in 0:n_segments-1)
    k_eff        = 1.0 / sum_inv_k
    α_avg        = u[1] / n_segments
    τ_transmitted = k_eff * sin(α_avg)

    k_bottom = spring_coeff * r_bottom
    α_bottom = abs(u[1]) * k_eff / k_bottom
    if abs(α_avg) >= 0.95 * π || α_bottom >= 0.95 * π
        τ_transmitted = 0.0
    end

    ω_ground = sqrt(max(τ_transmitted, 0.0) / p.k_mppt)
    return τ_transmitted * ω_ground
end
```

**Step 4: Add `trpt_ode_limited!`**

Append after `instantaneous_power` (3-arg variant):

```julia
"""
    trpt_ode_limited!(du, u, p::SystemParams, t)

3-state TRPT ODE with proportional elevation-angle power limiter.

State vector `u`:
- `u[1]` = α_tot : total TRPT twist angle (rad)
- `u[2]` = ω     : airborne rotor angular velocity (rad/s)
- `u[3]` = β     : elevation angle (rad) — dynamic state driven by power error

Control law:
    dβ/dt = clamp(kp_elev × (P_gen − p_rated_w), lower_rate, upper_rate)

where lower_rate = 0 when β ≤ β_min, upper_rate = 0 when β ≥ β_max.

Physics identical to `trpt_ode!` except `u[3]` replaces `p.elevation_angle` everywhere.
All existing 2-state ODEs are untouched.
"""
function trpt_ode_limited!(du, u, p::SystemParams, t)
    β = u[3]

    # Step 1 — Wind speed at hub altitude (using dynamic β)
    h     = hub_altitude(p.tether_length, β)
    v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h)

    # Steps 2–7 identical to trpt_ode! except β replaces p.elevation_angle ──────
    n_segments_inertia = p.n_rings + 1
    r_top_inertia  = p.trpt_hub_radius
    r_bot_inertia  = 2.0 * p.tether_length * p.trpt_rL_ratio / n_segments_inertia - r_top_inertia
    I_rings = sum(p.m_ring * (r_bot_inertia + i / p.n_rings * (r_top_inertia - r_bot_inertia))^2
                  for i in 1:p.n_rings)
    I_total = p.n_blades * p.m_blade * p.rotor_radius^2 + I_rings + p.i_pto

    n_segments = p.n_rings + 1
    r_top      = p.trpt_hub_radius
    r_bottom   = 2.0 * p.tether_length * p.trpt_rL_ratio / n_segments - r_top
    spring_coeff = p.e_modulus * π * (p.tether_diameter / 2)^2 * p.n_lines * p.trpt_rL_ratio
    sum_inv_k    = sum(1.0 / (spring_coeff * (r_bottom + i / (n_segments - 1) * (r_top - r_bottom)))
                       for i in 0:n_segments-1)
    k_eff = 1.0 / sum_inv_k

    α_avg         = u[1] / n_segments
    τ_transmitted = k_eff * sin(α_avg)

    k_bottom = spring_coeff * r_bottom
    α_bottom = abs(u[1]) * k_eff / k_bottom
    if abs(α_avg) >= 0.95 * π || α_bottom >= 0.95 * π
        τ_transmitted = 0.0
    end

    ω      = u[2]
    ω_safe = max(abs(ω), 0.1)
    λ_t    = ω * p.rotor_radius / max(v_hub, 0.1)
    P_aero = 0.5 * p.rho * v_hub^3 * π * p.rotor_radius^2 *
             cp_at_tsr(λ_t) * cos(β)^3                     # β from state, not p
    τ_aero = sign(ω) * P_aero / ω_safe

    V_a        = v_hub * (λ_t + sin(β))                    # β from state
    drag_force = 0.25 * 1.0 * p.tether_diameter * p.tether_length * p.rho * V_a^2
    τ_drag     = drag_force * p.rotor_radius * 0.5

    ω_ground = sqrt(max(τ_transmitted, 0.0) / p.k_mppt)

    du[1] = ω - ω_ground
    du[2] = (τ_aero - τ_drag - τ_transmitted) / I_total

    # Step 8 — Elevation-angle limiter (proportional control on ground power)
    P_gen    = τ_transmitted * ω_ground
    P_error  = P_gen - p.p_rated_w
    rate_raw = p.kp_elev * P_error
    lower_rate = (β <= p.β_min) ? 0.0 : -p.β_rate_max
    upper_rate = (β >= p.β_max) ? 0.0 :  p.β_rate_max
    du[3] = clamp(rate_raw, lower_rate, upper_rate)

    return nothing
end
```

**Step 5: Add `trpt_ode_wind_limited!`**

Append directly after `trpt_ode_limited!`:

```julia
"""
    trpt_ode_wind_limited!(du, u, pw, t)

Variable-wind 3-state variant. `pw = (p::SystemParams, wind_fn)`.
Physics and control law identical to `trpt_ode_limited!`.
"""
function trpt_ode_wind_limited!(du, u, pw, t)
    p, wind_fn = pw
    β = u[3]

    h       = hub_altitude(p.tether_length, β)
    v_ref_t = wind_fn(t)
    v_hub   = wind_at_altitude(v_ref_t, p.h_ref, h)

    n_segments_inertia = p.n_rings + 1
    r_top_inertia  = p.trpt_hub_radius
    r_bot_inertia  = 2.0 * p.tether_length * p.trpt_rL_ratio / n_segments_inertia - r_top_inertia
    I_rings = sum(p.m_ring * (r_bot_inertia + i / p.n_rings * (r_top_inertia - r_bot_inertia))^2
                  for i in 1:p.n_rings)
    I_total = p.n_blades * p.m_blade * p.rotor_radius^2 + I_rings + p.i_pto

    n_segments   = p.n_rings + 1
    r_top        = p.trpt_hub_radius
    r_bottom     = 2.0 * p.tether_length * p.trpt_rL_ratio / n_segments - r_top
    spring_coeff = p.e_modulus * π * (p.tether_diameter / 2)^2 * p.n_lines * p.trpt_rL_ratio
    sum_inv_k    = sum(1.0 / (spring_coeff * (r_bottom + i / (n_segments - 1) * (r_top - r_bottom)))
                       for i in 0:n_segments-1)
    k_eff = 1.0 / sum_inv_k

    α_avg         = u[1] / n_segments
    τ_transmitted = k_eff * sin(α_avg)

    k_bottom = spring_coeff * r_bottom
    α_bottom = abs(u[1]) * k_eff / k_bottom
    if abs(α_avg) >= 0.95 * π || α_bottom >= 0.95 * π
        τ_transmitted = 0.0
    end

    ω      = u[2]
    ω_safe = max(abs(ω), 0.1)
    λ_t    = ω * p.rotor_radius / max(v_hub, 0.1)
    P_aero = 0.5 * p.rho * v_hub^3 * π * p.rotor_radius^2 *
             cp_at_tsr(λ_t) * cos(β)^3
    τ_aero = sign(ω) * P_aero / ω_safe

    V_a        = v_hub * (λ_t + sin(β))
    drag_force = 0.25 * 1.0 * p.tether_diameter * p.tether_length * p.rho * V_a^2
    τ_drag     = drag_force * p.rotor_radius * 0.5

    ω_ground = sqrt(max(τ_transmitted, 0.0) / p.k_mppt)

    du[1] = ω - ω_ground
    du[2] = (τ_aero - τ_drag - τ_transmitted) / I_total

    P_gen    = τ_transmitted * ω_ground
    P_error  = P_gen - p.p_rated_w
    rate_raw = p.kp_elev * P_error
    lower_rate = (β <= p.β_min) ? 0.0 : -p.β_rate_max
    upper_rate = (β >= p.β_max) ? 0.0 :  p.β_rate_max
    du[3] = clamp(rate_raw, lower_rate, upper_rate)

    return nothing
end
```

**Step 6: Export new symbols from the module**

In `src/TRPTKiteTurbineSimulator.jl`, confirm `dynamics.jl` is included (it is). No additional exports needed — the module re-exports everything from included files.

**Step 7: Run full test suite**

```bash
julia --project=. test/runtests.jl
```

Expected: all tests pass including the new limiter tests. Count should be **384 + 5 = 389 tests**.

**Step 8: Commit**

```bash
git add src/dynamics.jl test/test_dynamics.jl
git commit -m "feat: add trpt_ode_limited!, trpt_ode_wind_limited!, instantaneous_power/3 (3-state elevation limiter)"
```

---

## Task 4: Add limiter trajectory support to build_trpt_scene

**Files:**
- Modify: `src/visualization.jl`

The `traj` named tuple gains an optional `beta` field. The scene must handle both:
- 2-state traj (no `beta` field) — existing behaviour
- 3-state traj (with `beta` field) — β displayed in HUD

**Step 1: Update `build_trpt_scene` to detect and normalise `beta`**

In `build_trpt_scene`, after `traj_norm = (t=..., v_hub=vhv)` (around line 453), add:

```julia
    # Normalise optional beta trace (3-state limiter trajectories include it)
    has_beta = hasproperty(traj, :beta)
    beta_vec = has_beta ? traj.beta : fill(rad2deg(p.elevation_angle), n_frames)
    traj_norm = (t=traj.t, alpha_tot=traj.alpha_tot, omega=traj.omega,
                 power_kw=traj.power_kw, v_hub=vhv, beta=beta_vec)
```

Add a `beta_max_obs` observable:

```julia
    beta_max_obs = Observable(maximum(beta_vec))
```

Pass `beta_max_obs` and `has_beta` through to `_build_hud!` and `_build_controls!`.

**Step 2: Add β row to `_build_hud!`**

In `_build_hud!`, find the power/ω readout rows. Add a β row after ω:

```julia
    # β row — only visible when trajectory has dynamic elevation
    β_lbl = Label(hud_layout[β_row, 1], "Elevation β"; halign=:left)
    β_val = Label(hud_layout[β_row, 2], "—°"; halign=:right)

    on(time_obs) do i
        β_deg = traj_obs[].beta[i]
        β_val.text[] = @sprintf("%.1f°", β_deg)
    end
```

**Step 3: Add "Auto limiter" toggle to `_build_controls!`**

After the existing `enable_toggle` row (layout[11]), add:

```julia
    # Auto limiter toggle (row 12-series additions)
    limiter_toggle = Toggle(layout[12, 1])
    Label(layout[12, 2], "Auto limiter (β control)"; halign=:left)
    Label(layout[13, 1], "Uses k_mppt + β_min below; β rises above rated power";
          halign=:left, fontsize=9, color=:grey60)
```

Shift subsequent layout rows down by 2 (rows 13+ become 15+).

**Step 4: Update re-run handler to branch on limiter toggle**

Replace the single `solve(ODEProblem(trpt_ode!, ...))` call with:

```julia
    if limiter_toggle.active[]
        # 3-state limiter mode: β evolves dynamically from β_min
        β0 = p_new.β_min
        sol_new = solve(ODEProblem(trpt_ode_limited!, [0.0, ω_warm, β0], (0.0, 120.0), p_new),
                        solver; kwargs...)
        n_frames_new = length(sol_new.t)
        v_hub_vec = [wind_at_altitude(p_new.v_wind_ref, p_new.h_ref,
                                      hub_altitude(p_new.tether_length, sol_new[3,i]))
                     for i in 1:n_frames_new]
        new_traj = (
            t         = sol_new.t,
            alpha_tot = sol_new[1, :],
            omega     = sol_new[2, :],
            power_kw  = [instantaneous_power(p_new, [sol_new[1,i], sol_new[2,i]], sol_new[3,i]) / 1000.0
                         for i in 1:n_frames_new],
            v_hub     = v_hub_vec,
            beta      = rad2deg.(sol_new[3, :]),
        )
    else
        # 2-state manual mode: fixed β from slider (existing behaviour)
        ω_warm = 4.1 * wind_at_altitude(v_new, p_new.h_ref,
                        hub_altitude(p_new.tether_length, β_new)) / p_new.rotor_radius
        sol_new = solve(ODEProblem(trpt_ode!, [0.0, ω_warm], (0.0, 120.0), p_new),
                        solver; kwargs...)
        n_frames_new = length(sol_new.t)
        v_hub_new = wind_at_altitude(p_new.v_wind_ref, p_new.h_ref,
                                     hub_altitude(p_new.tether_length, β_new))
        v_hub_vec = fill(v_hub_new, n_frames_new)
        new_traj = (
            t         = sol_new.t,
            alpha_tot = sol_new[1, :],
            omega     = sol_new[2, :],
            power_kw  = [instantaneous_power(p_new, [sol_new[1,i], sol_new[2,i]]) / 1000.0
                         for i in eachindex(sol_new.t)],
            v_hub     = v_hub_vec,
            beta      = fill(rad2deg(β_new), n_frames_new),
        )
    end
```

Note: also apply the warm-start `ω_warm` to the existing 2-state branch (replacing the hardcoded `1.0` initial condition — this is a long-overdue fix).

**Step 5: Run smoke test (headless)**

```bash
julia --project=. scripts/viz_smoke_test.jl
```

Expected: completes without error, saves PNG.

**Step 6: Commit**

```bash
git add src/visualization.jl
git commit -m "feat: add auto-limiter toggle and dynamic β HUD row to dashboard"
```

---

## Task 5: Create `scripts/run_wind_ramp_limited.jl`

**Files:**
- Create: `scripts/run_wind_ramp_limited.jl`

6 → 20 m/s ramp over 120 s. Warm-start at β_min, ω₀ for 6 m/s. The trajectory carries `beta` so the HUD shows elevation rising as rated power is reached.

```julia
# scripts/run_wind_ramp_limited.jl
# Wind ramp 6 → 20 m/s over 120 s with elevation-angle power limiter active.
# Demonstrates MPPT region (β flat at 23°) transitioning to limited region
# (β rising to shed excess power) as wind passes rated speed (~11 m/s at ~55 s).
#
# Usage:
#   julia --project=. scripts/run_wind_ramp_limited.jl

using DifferentialEquations, GLMakie, Statistics

include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/aerodynamics.jl")
include("../src/dynamics.jl")
include("../src/geometry.jl")
include("../src/force_analysis.jl")
include("../src/visualization.jl")

mkpath("output")

const SIM_DURATION = 120.0
const V_START      = 6.0     # m/s
const V_END        = 20.0    # m/s

p       = params_10kw()
wind_fn = wind_ramp(V_START, V_END, 0.0, SIM_DURATION)

# Warm-start: optimal TSR at the starting wind speed
h_hub = hub_altitude(p.tether_length, p.β_min)
v_hub0 = wind_at_altitude(V_START, p.h_ref, h_hub)
ω0 = 4.1 * v_hub0 / p.rotor_radius
u0 = [0.0, ω0, p.β_min]

println("Wind ramp: $(V_START) → $(V_END) m/s over $(SIM_DURATION) s")
println("Limiter: P_rated=$(p.p_rated_w/1000) kW, β_min=$(round(rad2deg(p.β_min),digits=1))°, β_max=$(round(rad2deg(p.β_max),digits=1))°")
println("Solving $(SIM_DURATION) s ODE...")

sol = solve(ODEProblem(trpt_ode_wind_limited!, u0, (0.0, SIM_DURATION), (p, wind_fn)),
            Tsit5(); reltol=1e-6, abstol=1e-6, saveat=1/30)
println("  $(length(sol.t)) frames  retcode: $(sol.retcode)")

v_hub_ts = [wind_at_altitude(wind_fn(t), p.h_ref,
             hub_altitude(p.tether_length, sol[3, i]))
            for (i, t) in enumerate(sol.t)]

traj = (
    t         = sol.t,
    alpha_tot = sol[1, :],
    omega     = sol[2, :],
    power_kw  = [instantaneous_power(p, [sol[1,i], sol[2,i]], sol[3,i]) / 1000.0
                 for i in eachindex(sol.t)],
    v_hub     = v_hub_ts,
    beta      = rad2deg.(sol[3, :]),
)

println("Peak power: $(round(maximum(traj.power_kw), digits=2)) kW")
println("β range:    $(round(minimum(traj.beta), digits=1))° – $(round(maximum(traj.beta), digits=1))°")
println("Mean power (last 30 s): $(round(mean(traj.power_kw[end-round(Int,30*30):end]), digits=2)) kW")

fig, _ = build_trpt_scene(p, traj)
display(fig)
readline()
```

**Step 1: Run the script (interactive test)**

```bash
julia --project=. scripts/run_wind_ramp_limited.jl
```

Expected:
- ODE solves successfully (retcode: Success)
- β trace rises from 23° and stabilises somewhere between 23° and 67°
- Peak power stays near P_rated (≤ ~12 kW with limiter settling time)
- Dashboard displays with β HUD row populated

**Step 2: Commit**

```bash
git add scripts/run_wind_ramp_limited.jl
git commit -m "feat: add run_wind_ramp_limited.jl (6→20 m/s, elevation power limiter)"
```

---

## Task 6: Update test suite — add limiter parameter tests

**Files:**
- Modify: `test/test_parameters.jl`

**Step 1: Add tests for new fields**

Append to `test/test_parameters.jl`:

```julia
@testset "SystemParams — limiter fields (10 kW preset)" begin
    p = params_10kw()
    @test p.p_rated_w  ≈ 10_000.0
    @test p.β_min      ≈ deg2rad(23.0) atol=0.001
    @test p.β_max      ≈ deg2rad(67.0) atol=0.001
    @test p.β_rate_max ≈ deg2rad(1.0)  atol=0.001
    @test p.kp_elev    ≈ 5e-5          atol=1e-6
    # β_min must be below β_max
    @test p.β_min < p.β_max
end

@testset "mass_scale — limiter fields scale correctly" begin
    p10 = params_10kw()
    p50 = mass_scale(p10, 10.0, 50.0)

    # Rated power scales linearly
    @test p50.p_rated_w ≈ p10.p_rated_w * 5.0

    # Angles do not scale
    @test p50.β_min      == p10.β_min
    @test p50.β_max      == p10.β_max
    @test p50.β_rate_max == p10.β_rate_max

    # kp_elev scales as 1/power_ratio (larger turbine: same °/s per fractional overpower)
    @test p50.kp_elev ≈ p10.kp_elev / 5.0 atol=1e-8
end
```

**Step 2: Run full test suite**

```bash
julia --project=. test/runtests.jl
```

Expected: **389 + 2 = 391 tests pass** (or more if additional assertions are counted separately).

**Step 3: Commit**

```bash
git add test/test_parameters.jl
git commit -m "test: add limiter field assertions to test_parameters.jl"
```

---

## Task 7: Update MEMORY.md

**Files:**
- Modify: `/home/rod/.claude/projects/-home-rod-Documents-GitHub-TRPTKiteTurbineJulia/memory/MEMORY.md`

Update the TODO section to mark TODO-1 complete. Update status line to reflect new test count and the new script. Add to the key files table:

```
| `scripts/run_wind_ramp_limited.jl` | 6→20 m/s ramp with elevation power limiter |
```

Add to function name reference:

```julia
# 3-state limiter ODEs
trpt_ode_limited!(du, u, p::SystemParams, t)             # u[3]=β, steady wind
trpt_ode_wind_limited!(du, u, pw, t)                      # u[3]=β, variable wind
instantaneous_power(p, u, β::Float64)                     # 3-state variant
```

---

## Final Verification

```bash
# Full test suite
julia --project=. test/runtests.jl

# Smoke test (headless)
julia --project=. scripts/viz_smoke_test.jl

# Power curve (verify no regression)
julia --project=. scripts/power_curve.jl
```

All 391 tests pass, smoke test saves PNG, power curve prints clean table.
