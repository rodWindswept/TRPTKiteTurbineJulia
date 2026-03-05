# TRPT 3D Visualization Engine Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend the Julia/GLMakie simulator with a high-fidelity 3D visualization engine:
inclined shaft geometry, rectangular rotor blades, per-element force coloring with auto-calibrated
FoS colourbar, elevation/azimuth/c_pto sliders, and on-demand ODE re-run with solver selection.

**Architecture:** Layered module refactor — `src/geometry.jl` (3D geometry), `src/force_analysis.jl`
(per-element forces), refactored `src/visualization.jl` (rendering + UI only). Existing `dynamics.jl`,
`parameters.jl`, `wind_profile.jl` untouched. TDD throughout: write failing test, implement, pass, commit.

**Tech Stack:** Julia 1.12, GLMakie, DifferentialEquations.jl, existing `SystemParams` struct.

---

## Background reading

Before starting, read these files to understand the codebase:
- `src/parameters.jl` — `SystemParams` struct and `params_10kw()` / `params_50kw()`
- `src/visualization.jl` — existing `compute_trpt_geometry()` and `build_trpt_scene()`
- `test/test_geometry.jl` — existing geometry tests (some will need updating)
- `docs/plans/2026-03-05-trpt-3d-viz-design.md` — approved design document

Key facts:
- `n_rings = 14`, `n_lines = 5`, `n_blades = 3` for 10 kW
- `trpt_hub_radius = 2.0 m` (rotor-end ring radius), `rotor_radius = 5.0 m` (blade tips)
- `tether_length = 30.0 m`, `elevation_angle = π/6` (30°)
- `r_bottom ≈ 0.96 m` (ground-end ring radius, computed from taper)
- Tests run with: `julia --project=. test/runtests.jl`

---

## Task 1: Create `src/geometry.jl` — inclined shaft geometry

**Files:**
- Create: `src/geometry.jl`
- Modify: `test/test_geometry.jl` (fix include path + update radius test)
- Modify: `src/visualization.jl` (remove `compute_trpt_geometry`, add `include("geometry.jl")`)

### Step 1: Write the failing tests in `test/test_geometry.jl`

Replace the file content with the following (the old radius test using `sqrt(x²+y²)` breaks
with inclined shaft — replace it with distance-from-ring-centre):

```julia
using Test

include("../src/parameters.jl")
include("../src/geometry.jl")   # was visualization.jl

@testset "TRPT geometry — node positions (inclined shaft)" begin
    p = params_10kw()
    alpha_tot = 0.3

    # Default shaft_dir derived from p.elevation_angle: [cos(β), 0, sin(β)]
    nodes = compute_trpt_geometry(p, alpha_tot)

    # Shape unchanged
    @test size(nodes) == (p.n_rings + 2, p.n_lines, 3)

    # Ground level at z = 0
    @test all(nodes[1, :, 3] .≈ 0.0)

    # Top level at hub altitude
    h_hub = p.tether_length * sin(p.elevation_angle)
    @test nodes[end, 1, 3] ≈ h_hub atol=1e-6

    # Top level x-position = tether_length * cos(elevation_angle)  [shaft leans downwind]
    x_hub = p.tether_length * cos(p.elevation_angle)
    @test nodes[end, 1, 1] ≈ x_hub atol=1e-2   # within 1 cm (ring centre offset, not node)

    # Each node is r_i away from its ring centre
    n_levels = p.n_rings + 2
    n_seg    = p.n_rings + 1
    l_seg    = p.tether_length / n_seg
    r_bottom = 2.0 * p.tether_length * p.trpt_rL_ratio / n_seg - p.trpt_hub_radius
    shaft_dir = [cos(p.elevation_angle), 0.0, sin(p.elevation_angle)]

    for i in 1:n_levels
        level_idx  = i - 1
        r_expected = r_bottom + level_idx / (n_levels - 1) * (p.trpt_hub_radius - r_bottom)
        ring_centre = level_idx * l_seg .* shaft_dir
        for j in 1:p.n_lines
            node = nodes[i, j, :]
            dist = sqrt(sum((node .- ring_centre).^2))
            @test dist ≈ r_expected atol=1e-6
        end
    end
end

@testset "Zero twist — lines radially aligned" begin
    p     = params_10kw()
    nodes = compute_trpt_geometry(p, 0.0)

    # With zero twist, angular position (in perp plane) is the same at ground and top
    # Use the y-coordinate of perp plane as proxy: with default shaft_dir in xz-plane,
    # the y-axis is one of the perp basis vectors.
    phi_ground = atan.(nodes[1, :, 2], nodes[1, :, 1])
    phi_top    = atan.(nodes[end, :, 2], nodes[end, :, 1])
    @test phi_ground ≈ phi_top atol=1e-6
end

@testset "Custom shaft_dir — vertical shaft" begin
    p         = params_10kw()
    shaft_dir = [0.0, 0.0, 1.0]
    nodes     = compute_trpt_geometry(p, 0.0, shaft_dir)

    # With vertical shaft, z increases linearly, x=y=0 at ring centres
    n_seg  = p.n_rings + 1
    l_seg  = p.tether_length / n_seg

    # Ground level z = 0
    @test all(nodes[1, :, 3] .≈ 0.0)

    # Top level z = tether_length (all altitude, no horizontal lean)
    @test nodes[end, 1, 3] ≈ p.tether_length atol=1e-6

    # Radii from shaft axis (x=0, y=0) equal r_i  [vertical shaft only — sqrt(x²+y²) valid]
    n_levels = p.n_rings + 2
    r_bottom = 2.0 * p.tether_length * p.trpt_rL_ratio / n_seg - p.trpt_hub_radius
    for i in 1:n_levels
        level_idx  = i - 1
        r_expected = r_bottom + level_idx / (n_levels - 1) * (p.trpt_hub_radius - r_bottom)
        for j in 1:p.n_lines
            r_actual = sqrt(nodes[i, j, 1]^2 + nodes[i, j, 2]^2)
            @test r_actual ≈ r_expected atol=1e-6
        end
    end
end
```

### Step 2: Run tests — confirm they fail (geometry.jl doesn't exist yet)

```bash
julia --project=. test/runtests.jl 2>&1 | tail -20
```
Expected: error loading `geometry.jl` (file not found).

### Step 3: Create `src/geometry.jl` with fixed `compute_trpt_geometry`

```julia
# src/geometry.jl
# 3D geometry for the TRPT tensegrity structure.
# Coordinate system: +x downwind, +y crosswind-right, +z vertical up.
# TRPT ground anchor at world origin (0,0,0).

const BLADE_INNER_FRAC = 0.30   # fraction of blade span inboard of the TRPT ring

"""
    compute_trpt_geometry(p, alpha_tot[, shaft_dir])

Compute 3D node positions for the TRPT tensegrity at a given total twist angle.

`shaft_dir` is a unit vector along the TRPT shaft axis (default: derived from
`p.elevation_angle` pointing downwind — `[cos(β), 0, sin(β)]`).

Returns array of size `(n_levels, n_lines, 3)`:
- Axis 1: ring levels, index 1 = ground anchor, index end = rotor end
- Axis 2: tether lines
- Axis 3: [x, y, z] world coordinates (metres)

Ring centres step along `shaft_dir` by `l_seg` per level. Each ring node is
offset radially in the plane perpendicular to `shaft_dir` at radius `r_i`.
"""
function compute_trpt_geometry(p::SystemParams, alpha_tot::Float64,
                                shaft_dir::Vector{Float64} = [cos(p.elevation_angle),
                                                              0.0,
                                                              sin(p.elevation_angle)])
    n_levels  = p.n_rings + 2
    n_seg     = p.n_rings + 1
    l_seg     = p.tether_length / n_seg
    alpha_seg = alpha_tot / n_seg

    r_bottom = 2.0 * p.tether_length * p.trpt_rL_ratio / n_seg - p.trpt_hub_radius

    # Perpendicular basis in the ring plane — use [0,0,1] as reference unless shaft is
    # nearly vertical, then use [0,1,0].
    ref  = abs(shaft_dir[3]) < 0.99 ? [0.0, 0.0, 1.0] : [0.0, 1.0, 0.0]
    perp1 = normalize(cross(shaft_dir, ref))   # first radial axis in ring plane
    perp2 = cross(shaft_dir, perp1)             # second radial axis (unit, perp to both)

    nodes = zeros(Float64, n_levels, p.n_lines, 3)

    for i in 1:n_levels
        level_idx   = i - 1
        r_i         = r_bottom + level_idx / (n_levels - 1) * (p.trpt_hub_radius - r_bottom)
        ring_centre = level_idx * l_seg .* shaft_dir
        theta_i     = level_idx * alpha_seg    # cumulative twist at this level

        for j in 1:p.n_lines
            phi  = theta_i + (j - 1) * (2π / p.n_lines)
            node = ring_centre .+ r_i .* (cos(phi) .* perp1 .+ sin(phi) .* perp2)
            nodes[i, j, :] = node
        end
    end

    return nodes
end

"""
    compute_blade_geometry(p, alpha_tot[, shaft_dir])

Compute 3D blade quad corners for all rotor blades at a given total twist angle.

Returns array of size `(n_blades, 4, 3)` — four corners per blade (quad outline).

Blades are locked in the rotor plane (perpendicular to `shaft_dir`) at the rotor end
of the TRPT. Each blade spans from `blade_inner_r` to `rotor_radius`, centred on the
shaft axis, with chord = 15% of `rotor_radius`.

The TRPT ring sits at `trpt_hub_radius`, which is approximately (1 - BLADE_INNER_FRAC)
of the total blade span from the outer tip — i.e., about 30% of the blade is inboard
of the ring.
"""
function compute_blade_geometry(p::SystemParams, alpha_tot::Float64,
                                 shaft_dir::Vector{Float64} = [cos(p.elevation_angle),
                                                               0.0,
                                                               sin(p.elevation_angle)])
    blade_span    = (p.rotor_radius - p.trpt_hub_radius) / (1.0 - BLADE_INNER_FRAC)
    blade_inner_r = p.trpt_hub_radius - BLADE_INNER_FRAC * blade_span
    blade_outer_r = p.rotor_radius
    chord         = p.rotor_radius * 0.15

    # Rotor centre (top end of TRPT shaft)
    n_seg       = p.n_rings + 1
    l_seg       = p.tether_length / n_seg
    rotor_centre = (p.n_rings + 1) * l_seg .* shaft_dir

    # Perpendicular basis (same logic as compute_trpt_geometry)
    ref   = abs(shaft_dir[3]) < 0.99 ? [0.0, 0.0, 1.0] : [0.0, 1.0, 0.0]
    perp1 = normalize(cross(shaft_dir, ref))
    perp2 = cross(shaft_dir, perp1)

    blades = zeros(Float64, p.n_blades, 4, 3)

    for b in 1:p.n_blades
        # Blade angle: cumulative twist at rotor end + even angular spacing
        phi_b = alpha_tot + (b - 1) * (2π / p.n_blades)

        # Blade radial direction in ring plane
        blade_dir = cos(phi_b) .* perp1 .+ sin(phi_b) .* perp2

        # Chord direction — perpendicular to blade_dir and shaft_dir
        chord_dir = cross(shaft_dir, blade_dir)

        inner_centre = rotor_centre .+ blade_inner_r .* blade_dir
        outer_centre = rotor_centre .+ blade_outer_r .* blade_dir

        blades[b, 1, :] = inner_centre .+ (chord / 2) .* chord_dir
        blades[b, 2, :] = inner_centre .- (chord / 2) .* chord_dir
        blades[b, 3, :] = outer_centre .- (chord / 2) .* chord_dir
        blades[b, 4, :] = outer_centre .+ (chord / 2) .* chord_dir
    end

    return blades
end

"""
    world_ground_plane(half_extent=40.0, step=5.0)

Return a vector of line segments forming a grid at z = 0.
Each element is a tuple `(xs, ys, zs)` of two-point coordinate vectors.
"""
function world_ground_plane(half_extent::Float64 = 40.0, step::Float64 = 5.0)
    lines = Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}[]
    coords = -half_extent:step:half_extent
    for c in coords
        push!(lines, ([-half_extent, half_extent], [c, c],              [0.0, 0.0]))
        push!(lines, ([c, c],              [-half_extent, half_extent], [0.0, 0.0]))
    end
    return lines
end
```

Note: `normalize` and `cross` come from `LinearAlgebra`, which GLMakie re-exports.
If running geometry.jl standalone in tests, add `using LinearAlgebra` at the top.

### Step 4: Remove `compute_trpt_geometry` from `src/visualization.jl`

Delete lines 7–50 of `visualization.jl` (the `compute_trpt_geometry` function and its
docstring). Add `include("geometry.jl")` after `using GLMakie` at the top.

### Step 5: Update `src/TRPTKiteTurbineSimulator.jl`

Add `include("geometry.jl")` before `include("visualization.jl")`:

```julia
module TRPTKiteTurbineSimulator

include("parameters.jl")
include("wind_profile.jl")
include("dynamics.jl")
include("geometry.jl")
include("visualization.jl")

end
```

### Step 6: Run full test suite — all tests pass

```bash
julia --project=. test/runtests.jl
```
Expected: all existing tests pass + 3 new geometry testsets pass.

### Step 7: Commit

```bash
git add src/geometry.jl src/visualization.jl src/TRPTKiteTurbineSimulator.jl test/test_geometry.jl
git commit -m "feat: extract geometry.jl with inclined shaft fix and blade geometry"
```

---

## Task 2: Create `src/force_analysis.jl`

**Files:**
- Create: `src/force_analysis.jl`
- Create: `test/test_force_analysis.jl`
- Modify: `test/runtests.jl` (add include)
- Modify: `src/TRPTKiteTurbineSimulator.jl` (add include)

### Step 1: Write failing tests in `test/test_force_analysis.jl`

```julia
using Test

include("../src/parameters.jl")
include("../src/wind_profile.jl")
include("../src/dynamics.jl")
include("../src/force_analysis.jl")

@testset "ForceState — structure" begin
    p     = params_10kw()
    n_seg = p.n_rings + 1
    h_hub = hub_altitude(p.tether_length, p.elevation_angle)
    v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h_hub)
    u     = [0.5, 1.8]   # typical [alpha_tot, omega] at rated conditions

    fs = element_forces(p, u, v_hub)

    @test fs isa ForceState
    @test length(fs.tether_tension)   == n_seg
    @test length(fs.ring_compression) == p.n_rings
    @test fs.tau_transmitted          >= 0.0
    @test fs.tau_aero                 >  0.0
    @test fs.tau_drag                 >= 0.0
end

@testset "ForceState — physical bounds" begin
    p     = params_10kw()
    h_hub = hub_altitude(p.tether_length, p.elevation_angle)
    v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h_hub)
    u     = [0.5, 1.8]

    fs = element_forces(p, u, v_hub)

    # All tensions and compressions must be non-negative
    @test all(fs.tether_tension   .>= 0.0)
    @test all(fs.ring_compression .>= 0.0)

    # Tether SWL for 3mm Dyneema is ~3500 N; in normal operation tensions << SWL
    @test all(fs.tether_tension .< 3500.0)

    # Ring compressions must be positive and below 500 N (CFRP conservative estimate)
    @test all(fs.ring_compression .< 500.0)
end

@testset "ForceState — zero omega gives zero centrifugal" begin
    p     = params_10kw()
    h_hub = hub_altitude(p.tether_length, p.elevation_angle)
    v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h_hub)
    u0    = [0.0, 0.0]   # no twist, no rotation

    fs0 = element_forces(p, u0, v_hub)

    # With zero twist, transmitted torque = k_eff * sin(0) = 0 → no torsional tension
    # Remaining tension is purely gravitational (m_above * g * sin(β))
    @test fs0.tau_transmitted ≈ 0.0 atol=1e-6
    @test all(fs0.tether_tension .>= 0.0)
end

@testset "run_force_scan — returns max over trajectory" begin
    p     = params_10kw()
    h_hub = hub_altitude(p.tether_length, p.elevation_angle)
    v_hub = wind_at_altitude(p.v_wind_ref, p.h_ref, h_hub)

    # Fake 3-frame trajectory
    u_frames   = [[0.1, 1.0], [0.5, 1.8], [0.3, 1.5]]
    v_hub_vec  = fill(v_hub, 3)

    T_max, C_max = run_force_scan(p, u_frames, v_hub_vec)

    @test T_max > 0.0
    @test C_max > 0.0
    # Max must be >= any individual frame
    for (u, v) in zip(u_frames, v_hub_vec)
        fs = element_forces(p, u, v)
        @test T_max >= maximum(fs.tether_tension)
        @test C_max >= maximum(fs.ring_compression)
    end
end
```

### Step 2: Run tests — confirm they fail

```bash
julia --project=. test/runtests.jl 2>&1 | grep -A3 "force_analysis"
```
Expected: error — `ForceState` not defined.

### Step 3: Create `src/force_analysis.jl`

```julia
# src/force_analysis.jl
# Per-element structural force analysis derived analytically from ODE state.
# All forces in SI units (N, N·m).

const TETHER_SWL = 3500.0   # Dyneema 3mm safe working load (N) — DRR §5.2
const RING_SWL   = 500.0    # CFRP tube ring strut conservative buckling limit (N)

"""
    ForceState

Per-element structural forces at one ODE state.
"""
struct ForceState
    tether_tension   :: Vector{Float64}   # axial tension per TRPT segment (N), length n_seg
    ring_compression :: Vector{Float64}   # hoop compression per ring (N), length n_rings
    tau_transmitted  :: Float64           # torque through TRPT (N·m)
    tau_aero         :: Float64           # aerodynamic torque (N·m)
    tau_drag         :: Float64           # tether drag torque (N·m)
end

"""
    element_forces(p, u, v_hub) -> ForceState

Derive per-element structural forces analytically from ODE state `u = [alpha_tot, omega]`
and hub wind speed `v_hub` (m/s).

Tether tension per segment i:
    T_i = τ_torsion_i + τ_centrifugal_i + τ_gravity_i

where:
- τ_torsion_i    = (k_i × sin(α_i)) / r_i    (torsional torque arm)
- τ_centrifugal_i = m_seg × ω² × r_i          (centrifugal load on tether segment)
- τ_gravity_i     = m_above_i × g × sin(β)    (gravitational component along shaft)

Ring hoop compression from polygon geometry:
    C_i = T_i × r_i / (2 × sin(π / n_lines))
"""
function element_forces(p::SystemParams, u::Vector{Float64}, v_hub::Float64)::ForceState
    n_seg    = p.n_rings + 1
    r_top    = p.trpt_hub_radius
    r_bottom = 2.0 * p.tether_length * p.trpt_rL_ratio / n_seg - r_top
    l_seg    = p.tether_length / n_seg

    spring_coeff = p.e_modulus * π * (p.tether_diameter / 2)^2 * p.n_lines * p.trpt_rL_ratio

    alpha_tot = u[1]
    omega     = u[2]
    g         = 9.81

    # Mass of one tether segment
    tether_line_mass_per_m = p.rho * π * (p.tether_diameter / 2)^2  # linear density (kg/m)
    m_seg_tether = tether_line_mass_per_m * l_seg * p.n_lines         # all lines in segment

    # Cumulative twist per segment (uniform distribution)
    alpha_per_seg = alpha_tot / n_seg

    tether_tension   = zeros(Float64, n_seg)
    ring_compression = zeros(Float64, p.n_rings)

    # τ_transmitted (same calculation as dynamics.jl step 4)
    sum_inv_k = sum(1.0 / (spring_coeff * (r_bottom + i / (n_seg - 1) * (r_top - r_bottom)))
                    for i in 0:n_seg-1)
    k_eff = 1.0 / sum_inv_k
    alpha_avg = alpha_tot / n_seg
    tau_transmitted = k_eff * sin(alpha_avg)

    # Aerodynamic torque
    omega_safe = max(abs(omega), 0.1)
    P_aero     = 0.5 * p.rho * v_hub^3 * π * p.rotor_radius^2 * p.cp * cos(p.elevation_angle)^3
    tau_aero   = sign(omega) * P_aero / omega_safe

    # Tether drag torque
    lambda_t   = omega * p.rotor_radius / max(v_hub, 0.1)
    V_a        = v_hub * (lambda_t + sin(p.elevation_angle))
    drag_force = 0.25 * 1.0 * p.tether_diameter * p.tether_length * p.rho * V_a^2
    tau_drag   = drag_force * p.rotor_radius * 0.5

    for i in 0:n_seg-1
        # Radius at segment i (0 = ground, n_seg-1 = rotor)
        r_i = r_bottom + i / (n_seg - 1) * (r_top - r_bottom)
        k_i = spring_coeff * r_i

        # Torsional component: local torque / moment arm
        tau_torsion_i = k_i * sin(alpha_per_seg) / r_i

        # Centrifugal component: segment mass × ω² × r_i
        tau_centrifugal_i = m_seg_tether * omega^2 * r_i

        # Gravitational component: mass above this segment projected along shaft
        # mass_above = rings above + blades + remaining tether segments above
        n_above      = n_seg - 1 - i
        m_rings_above = n_above * p.m_ring
        m_blades_above = (i == n_seg - 1) ? p.n_blades * p.m_blade : 0.0
        m_tether_above = n_above * m_seg_tether
        m_above       = m_rings_above + m_blades_above + m_tether_above
        tau_gravity_i = m_above * g * sin(p.elevation_angle)

        tether_tension[i + 1] = max(0.0, tau_torsion_i + tau_centrifugal_i + tau_gravity_i)

        # Ring hoop compression from polygon tension
        if i < p.n_rings
            ring_compression[i + 1] = tether_tension[i + 1] * r_i /
                                       (2.0 * sin(π / p.n_lines))
        end
    end

    return ForceState(tether_tension, ring_compression, tau_transmitted, tau_aero, tau_drag)
end

"""
    run_force_scan(p, u_frames, v_hub_vec) -> (T_max, C_max)

Scan a full pre-computed trajectory to find the maximum per-element forces.
`u_frames` is a vector of state vectors `[alpha_tot, omega]`.
`v_hub_vec` is a vector of hub wind speeds (m/s), one per frame.

Returns `(T_max_run, C_max_run)` — the peak tether tension and ring compression
observed across all frames, used to calibrate the force colourbar.
"""
function run_force_scan(p::SystemParams,
                        u_frames::Vector{Vector{Float64}},
                        v_hub_vec::Vector{Float64})
    T_max = 0.0
    C_max = 0.0
    for (u, v) in zip(u_frames, v_hub_vec)
        fs    = element_forces(p, u, v)
        T_max = max(T_max, maximum(fs.tether_tension))
        C_max = max(C_max, maximum(fs.ring_compression))
    end
    return T_max, C_max
end
```

### Step 4: Update `test/runtests.jl`

Add `include("test_force_analysis.jl")` inside the `@testset` block.

### Step 5: Update `src/TRPTKiteTurbineSimulator.jl`

Add `include("force_analysis.jl")` after `include("geometry.jl")`.

### Step 6: Run full test suite

```bash
julia --project=. test/runtests.jl
```
Expected: all tests pass including 4 new force_analysis testsets.

### Step 7: Commit

```bash
git add src/force_analysis.jl test/test_force_analysis.jl test/runtests.jl src/TRPTKiteTurbineSimulator.jl
git commit -m "feat: add force_analysis.jl with per-element tension and compression"
```

---

## Task 3: Refactor `src/visualization.jl` — geometry observables and scene upgrade

This is the largest task. Work through it section by section.

**Files:**
- Modify: `src/visualization.jl` (full rewrite of `build_trpt_scene`)
- Modify: `scripts/run_dashboard_10kw.jl` (add geometry.jl include)
- Modify: `scripts/run_dashboard_50kw.jl` (add geometry.jl include)

**Note:** There are no automated tests for the visual scene (GLMakie requires a display).
Smoke-test manually by running `julia --project=. scripts/run_dashboard_10kw.jl` and
verifying the window opens with the new layout.

### Step 1: Update dashboard script includes

In `scripts/run_dashboard_10kw.jl` and `scripts/run_dashboard_50kw.jl`, add:
```julia
include("../src/geometry.jl")
include("../src/force_analysis.jl")
```
after the existing `include("../src/dynamics.jl")` line.

### Step 2: Rewrite `src/visualization.jl`

Replace the entire file with the following:

```julia
# src/visualization.jl
# GLMakie 3D visualization of the TRPT tensegrity structure.
# Depends on: geometry.jl, force_analysis.jl, parameters.jl, dynamics.jl
# Strategy: pre-compute ODE trajectory + force scan, then animate at 30 FPS.

using GLMakie
using LinearAlgebra
using DifferentialEquations

# ── Colour helpers ─────────────────────────────────────────────────────────────

"""Map a scalar `v` in [0, v_max] to a blue→red RGB colour."""
function _force_color(v::Float64, v_max::Float64)
    t = v_max > 0.0 ? clamp(v / v_max, 0.0, 1.0) : 0.0
    return RGBf(t, 0.0, 1.0 - t)   # blue (t=0) → red (t=1)
end

# ── 3D axis builder ────────────────────────────────────────────────────────────

function _build_3d_axes!(fig, position, p, traj, shaft_dir_obs, T_max_obs, C_max_obs,
                          time_obs, v_hub_vec)
    ax3d = Axis3(fig[position...],
                 title  = "TRPT Kite Turbine — Live Geometry",
                 xlabel = "X downwind (m)", ylabel = "Y crosswind (m)", zlabel = "Altitude (m)",
                 aspect = :data)

    # Ground plane grid (static)
    for (xs, ys, zs) in world_ground_plane()
        lines!(ax3d, xs, ys, zs; color=(:grey, 0.3), linewidth=0.5)
    end

    # Wind direction arrow (static, updates with azimuth slider via shaft_dir_obs)
    wind_arrow = @lift begin
        sd = $shaft_dir_obs
        # Arrow in ground plane, pointing in wind direction
        ([0.0, sd[1] * 15.0], [0.0, sd[2] * 15.0], [0.5, 0.5])
    end
    lines!(ax3d, @lift($wind_arrow[1]), @lift($wind_arrow[2]), @lift($wind_arrow[3]);
           color=:darkorange, linewidth=3.0)

    # Nodes observable (reacts to time AND shaft_dir)
    nodes_obs = @lift begin
        compute_trpt_geometry(p, traj.alpha_tot[$time_obs], $shaft_dir_obs)
    end

    # Blades observable
    blades_obs = @lift begin
        compute_blade_geometry(p, traj.alpha_tot[$time_obs], $shaft_dir_obs)
    end

    # Force colours per tether line (n_lines lines, each with n_seg segments coloured)
    # We draw each tether line as a series of coloured segments.
    n_seg    = p.n_rings + 1
    n_levels = p.n_rings + 2

    for j in 1:p.n_lines
        # Per-frame node positions for line j
        xs = @lift $nodes_obs[:, j, 1]
        ys = @lift $nodes_obs[:, j, 2]
        zs = @lift $nodes_obs[:, j, 3]

        # Per-frame colour for line j — colour changes per segment
        # GLMakie `lines!` with per-vertex colours: pass a vector of colours
        line_colors = @lift begin
            u     = [$traj.alpha_tot[$time_obs], traj.omega[$time_obs]]
            v_hub = v_hub_vec isa AbstractVector ? v_hub_vec[$time_obs] : v_hub_vec
            fs    = element_forces(p, u, v_hub)
            T_max = $T_max_obs
            # One colour per level (n_levels vertices) — repeat segment colour for both endpoints
            colors = Vector{RGBf}(undef, n_levels)
            colors[1] = _force_color(fs.tether_tension[1], T_max)
            for i in 2:n_levels
                seg_idx = min(i - 1, n_seg)
                colors[i] = _force_color(fs.tether_tension[seg_idx], T_max)
            end
            colors
        end

        lines!(ax3d, xs, ys, zs; color=line_colors, linewidth=1.5)
    end

    # Ring polygons — colour by ring compression
    for i in 2:(n_levels - 1)
        ring_x = @lift [$nodes_obs[i, :, 1]; $nodes_obs[i, 1, 1]]
        ring_y = @lift [$nodes_obs[i, :, 2]; $nodes_obs[i, 1, 2]]
        ring_z = @lift [$nodes_obs[i, :, 3]; $nodes_obs[i, 1, 3]]

        ring_color = @lift begin
            u     = [traj.alpha_tot[$time_obs], traj.omega[$time_obs]]
            v_hub = v_hub_vec isa AbstractVector ? v_hub_vec[$time_obs] : v_hub_vec
            fs    = element_forces(p, u, v_hub)
            C_max = $C_max_obs
            ring_idx = i - 1   # rings are indices 2..(n_levels-1)
            _force_color(fs.ring_compression[ring_idx], C_max)
        end

        lines!(ax3d, ring_x, ring_y, ring_z; color=ring_color, linewidth=1.0)
    end

    # Rotor ring (always red — structural landmark)
    rotor_x = @lift [$nodes_obs[end, :, 1]; $nodes_obs[end, 1, 1]]
    rotor_y = @lift [$nodes_obs[end, :, 2]; $nodes_obs[end, 1, 2]]
    rotor_z = @lift [$nodes_obs[end, :, 3]; $nodes_obs[end, 1, 3]]
    lines!(ax3d, rotor_x, rotor_y, rotor_z; color=:firebrick, linewidth=3.5)

    # Blades
    for b in 1:p.n_blades
        bx = @lift [$blades_obs[b, [1,2,3,4,1], 1]...]
        by = @lift [$blades_obs[b, [1,2,3,4,1], 2]...]
        bz = @lift [$blades_obs[b, [1,2,3,4,1], 3]...]
        lines!(ax3d, bx, by, bz; color=:steelblue, linewidth=2.0)
    end

    # Ground anchor
    scatter!(ax3d, [0.0], [0.0], [0.0]; color=:green, markersize=20)

    return ax3d
end

# ── Telemetry HUD builder ──────────────────────────────────────────────────────

function _build_hud!(parent_layout, p, traj, time_obs, T_max_obs, C_max_obs, v_hub_vec)
    Label(parent_layout[1, 1], "Telemetry"; fontsize=16, font=:bold, halign=:left)
    time_label   = Label(parent_layout[2, 1], "t = 0.00 s";        halign=:left)
    power_label  = Label(parent_layout[3, 1], "P = 0.00 kW";       halign=:left)
    omega_label  = Label(parent_layout[4, 1], "ω = 0.00 rad/s";    halign=:left)
    twist_label  = Label(parent_layout[5, 1], "α_seg = 0.00°";     halign=:left)
    margin_label = Label(parent_layout[6, 1], "Collapse margin: 100%"; halign=:left)
    wind_label   = Label(parent_layout[7, 1], "V_hub = 0.00 m/s";  halign=:left)

    Label(parent_layout[8, 1], ""; halign=:left)   # spacer

    # Force colourbar — tether tension
    Label(parent_layout[9, 1], "Tether tension (N)"; halign=:left, fontsize=12, font=:bold)
    tether_bar_label = Label(parent_layout[10, 1], "0 → 0 N  |  FoS: —"; halign=:left)

    # Force colourbar — ring compression
    Label(parent_layout[11, 1], "Ring compression (N)"; halign=:left, fontsize=12, font=:bold)
    ring_bar_label = Label(parent_layout[12, 1], "0 → 0 N  |  FoS: —"; halign=:left)

    n_seg = p.n_rings + 1

    on(time_obs) do frame_idx
        t       = traj.t[frame_idx]
        alpha   = traj.alpha_tot[frame_idx]
        omega   = traj.omega[frame_idx]
        power   = traj.power_kw[frame_idx]
        alpha_s = rad2deg(alpha / n_seg)
        margin  = max(0.0, (1.0 - abs(alpha / n_seg) / π) * 100.0)
        v_hub   = v_hub_vec isa AbstractVector ? v_hub_vec[frame_idx] : v_hub_vec

        time_label.text[]   = "t = $(round(t,      digits=2)) s"
        power_label.text[]  = "P = $(round(power,  digits=2)) kW"
        omega_label.text[]  = "ω = $(round(omega,  digits=3)) rad/s"
        twist_label.text[]  = "α_seg = $(round(alpha_s, digits=1))°"
        margin_label.text[] = "Collapse margin: $(round(margin, digits=1))%"
        wind_label.text[]   = "V_hub = $(round(v_hub, digits=2)) m/s"

        T_max = T_max_obs[]
        C_max = C_max_obs[]
        T_fos = T_max > 0 ? round(TETHER_SWL / T_max, digits=1) : Inf
        C_fos = C_max > 0 ? round(RING_SWL   / C_max, digits=1) : Inf
        tether_bar_label.text[] = "0 → $(round(T_max, digits=0)) N  |  SWL=$(TETHER_SWL)N  FoS=$(T_fos)"
        ring_bar_label.text[]   = "0 → $(round(C_max, digits=0)) N  |  SWL=$(RING_SWL)N  FoS=$(C_fos)"
    end

    return (time_label, power_label, omega_label, twist_label, margin_label,
            wind_label, tether_bar_label, ring_bar_label)
end

# ── Controls builder ───────────────────────────────────────────────────────────

function _build_controls!(parent_layout, p, traj_ref, time_obs_ref,
                           shaft_dir_obs, T_max_obs, C_max_obs, n_frames)

    Label(parent_layout[1, 1], "Controls"; fontsize=14, font=:bold, halign=:left)

    # Elevation angle slider (0–75°)
    Label(parent_layout[2, 1], "Elevation β (°)"; halign=:left)
    elev_slider = Slider(parent_layout[3, 1],
                         range=0.0:1.0:75.0,
                         startvalue=rad2deg(p.elevation_angle))

    # Wind azimuth slider (0–360°)
    Label(parent_layout[4, 1], "Wind azimuth φ (°)"; halign=:left)
    azimuth_slider = Slider(parent_layout[5, 1],
                            range=0.0:1.0:360.0,
                            startvalue=0.0)

    # Update shaft_dir when either angle slider changes
    on(elev_slider.value) do _
        β  = deg2rad(elev_slider.value[])
        φw = deg2rad(azimuth_slider.value[])
        shaft_dir_obs[] = [cos(φw)*cos(β), sin(φw)*cos(β), sin(β)]
    end
    on(azimuth_slider.value) do _
        β  = deg2rad(elev_slider.value[])
        φw = deg2rad(azimuth_slider.value[])
        shaft_dir_obs[] = [cos(φw)*cos(β), sin(φw)*cos(β), sin(β)]
    end

    # Time slider
    Label(parent_layout[6, 1], "Time"; halign=:left)
    time_slider = Slider(parent_layout[7, 1], range=1:n_frames, startvalue=1)
    time_obs_ref[] = time_slider.value

    # Play / Pause
    play_row = GridLayout(parent_layout[8, 1])
    play_button = Button(play_row[1, 1], label="▶ Play")

    # Solver menu
    Label(play_row[1, 2], "Solver:"; halign=:right)
    solver_menu = Menu(play_row[1, 3],
                       options=["Tsit5", "RK4", "Euler"],
                       default="Tsit5")

    # Animation loop
    is_playing = Observable(false)
    on(play_button.clicks) do _
        is_playing[] = !is_playing[]
        play_button.label[] = is_playing[] ? "⏸ Pause" : "▶ Play"
    end
    @async while true
        if is_playing[]
            next_frame = min(time_slider.value[] + 1, n_frames)
            set_close_to!(time_slider, next_frame)
            if next_frame == n_frames
                is_playing[] = false
                play_button.label[] = "▶ Play"
            end
        end
        sleep(1 / 30)
    end

    # Dynamic torque panel
    Label(parent_layout[9, 1], "Dynamic Torque Mode"; fontsize=13, font=:bold, halign=:left)
    enable_check  = Toggle(parent_layout[10, 1])
    Label(parent_layout[10, 2], "Enable"; halign=:left)

    Label(parent_layout[11, 1], "c_pto (N·m·s/rad)"; halign=:left)
    c_pto_slider  = Slider(parent_layout[12, 1],
                           range=exp10.(range(log10(100.0), log10(50000.0), length=200)),
                           startvalue=p.c_pto)
    c_pto_label   = Label(parent_layout[13, 1],
                          "$(round(p.c_pto, digits=0)) N·m·s/rad"; halign=:left)
    on(c_pto_slider.value) do v
        c_pto_label.text[] = "$(round(v, digits=0)) N·m·s/rad"
    end

    rerun_button  = Button(parent_layout[14, 1], label="Re-run ODE")

    on(rerun_button.clicks) do _
        if !enable_check.active[]
            return
        end

        c_pto_new = c_pto_slider.value[]
        # Rebuild params with new c_pto (all other fields unchanged)
        p_new = SystemParams(p.rho, p.v_wind_ref, p.h_ref, p.elevation_angle,
                             p.rotor_radius, p.tether_length, p.trpt_hub_radius,
                             p.trpt_rL_ratio, p.n_lines, p.tether_diameter,
                             p.e_modulus, p.n_rings, p.m_ring, p.n_blades,
                             p.m_blade, p.cp, p.i_pto, c_pto_new)

        solver = if solver_menu.selection[] == "Tsit5"
            Tsit5()
        elseif solver_menu.selection[] == "RK4"
            RK4()
        else
            Euler()
        end

        dt_arg = solver_menu.selection[] == "Euler" ? (dt=0.01,) : NamedTuple()

        sol_new = solve(ODEProblem(trpt_ode!, [0.0, 1.0], (0.0, 120.0), p_new),
                        solver; reltol=1e-6, abstol=1e-6, saveat=1/30, dt_arg...)

        new_traj = (
            t         = sol_new.t,
            alpha_tot = sol_new[1, :],
            omega     = sol_new[2, :],
            power_kw  = [instantaneous_power(p_new, [sol_new[1,i], sol_new[2,i]]) / 1000.0
                         for i in eachindex(sol_new.t)],
            v_hub     = traj_ref[].v_hub,
        )

        # Re-scan forces for new trajectory
        v_hub_scalar = new_traj.v_hub isa AbstractVector ? new_traj.v_hub[1] : new_traj.v_hub
        u_frames  = [[new_traj.alpha_tot[i], new_traj.omega[i]] for i in eachindex(new_traj.t)]
        v_hub_vec_new = fill(v_hub_scalar, length(new_traj.t))
        T_max_new, C_max_new = run_force_scan(p_new, u_frames, v_hub_vec_new)

        # Update observables atomically
        traj_ref[]  = new_traj
        T_max_obs[] = T_max_new
        C_max_obs[] = C_max_new
        set_close_to!(time_slider, 1)
    end

    return (time_slider, elev_slider, azimuth_slider, c_pto_slider, solver_menu)
end

# ── Main scene builder ─────────────────────────────────────────────────────────

"""
    build_trpt_scene(p, traj)

Build the GLMakie Figure for TRPT visualization and interactive control.

`traj` fields: `.t`, `.alpha_tot`, `.omega`, `.power_kw`, `.v_hub` (scalar or vector).

Returns `(fig, time_obs)`.
"""
function build_trpt_scene(p::SystemParams, traj)
    n_frames  = length(traj.t)
    v_hub_vec = isa(traj.v_hub, AbstractVector) ? traj.v_hub : fill(traj.v_hub, n_frames)

    # Pre-scan forces for initial colourbar calibration
    u_frames = [[traj.alpha_tot[i], traj.omega[i]] for i in 1:n_frames]
    T_max_run, C_max_run = run_force_scan(p, u_frames, v_hub_vec)

    fig = Figure(size=(1600, 900))

    # Shared observables
    shaft_dir_obs = Observable([cos(p.elevation_angle), 0.0, sin(p.elevation_angle)])
    T_max_obs     = Observable(T_max_run)
    C_max_obs     = Observable(C_max_run)
    time_obs      = Observable(1)
    time_obs_ref  = Ref{Any}(time_obs)
    traj_ref      = Ref(traj)

    # Layout: 3D viewport (left) | HUD + controls (right, 320px)
    colsize!(fig.layout, 2, Fixed(320))

    right_panel = GridLayout(fig[1, 2])
    hud_layout  = GridLayout(right_panel[1, 1])
    ctrl_layout = GridLayout(right_panel[2, 1])

    _build_3d_axes!(fig, (1, 1), p, traj, shaft_dir_obs, T_max_obs, C_max_obs,
                    time_obs, v_hub_vec)
    _build_hud!(hud_layout, p, traj, time_obs, T_max_obs, C_max_obs, v_hub_vec)
    _build_controls!(ctrl_layout, p, traj_ref, time_obs_ref, shaft_dir_obs,
                     T_max_obs, C_max_obs, n_frames)

    return fig, time_obs
end
```

### Step 3: Update remaining dashboard scripts

All variable-wind scripts (`run_wind_ramp.jl`, `run_gust.jl`, `run_turbulent.jl`) also
`include("../src/visualization.jl")` — add the same geometry/force_analysis includes to each.

### Step 4: Smoke test — run the 10 kW dashboard

```bash
julia --project=. scripts/run_dashboard_10kw.jl
```

Expected: GLMakie window opens showing inclined shaft, coloured tether lines, blades,
ground plane grid, wind arrow, HUD with FoS labels, Controls panel with all sliders.
Verify:
- Elevation slider tilts the whole shaft
- Wind azimuth rotates the shaft and wind arrow in the ground plane
- Play button animates the shaft twist and blade rotation
- Blades stay in the rotor plane at all elevation angles

### Step 5: Commit

```bash
git add src/visualization.jl src/geometry.jl scripts/run_dashboard_10kw.jl \
        scripts/run_dashboard_50kw.jl scripts/run_wind_ramp.jl \
        scripts/run_gust.jl scripts/run_turbulent.jl
git commit -m "feat: refactor visualization with inclined shaft, blades, force coloring, dynamic torque panel"
```

---

## Task 4: Extend CSV export with force columns

**Files:**
- Modify: `src/main.jl`

### Step 1: Read `src/main.jl` and locate the CSV export block

The file writes a `DataFrame` and saves to `output/`. Find the `DataFrame(...)` call.

### Step 2: Extend the DataFrame with force columns

After computing `traj`, add a force scan and extend the CSV:

```julia
# Force scan
u_frames  = [[sol[1,i], sol[2,i]] for i in eachindex(sol.t)]
v_hub_all = fill(v_hub, length(sol.t))
T_max_run, C_max_run = run_force_scan(p, u_frames, v_hub_all)

force_states = [element_forces(p, u_frames[i], v_hub_all[i]) for i in eachindex(sol.t)]
n_seg = p.n_rings + 1

df = DataFrame(
    t          = sol.t,
    alpha_tot  = sol[1, :],
    omega      = sol[2, :],
    power_kw   = traj.power_kw,
    v_hub      = v_hub_all,
    tau_aero        = [fs.tau_aero        for fs in force_states],
    tau_drag        = [fs.tau_drag        for fs in force_states],
    tau_transmitted = [fs.tau_transmitted for fs in force_states],
    T_max_run  = fill(T_max_run,  length(sol.t)),
    C_max_run  = fill(C_max_run,  length(sol.t)),
)

# Add per-segment tether tension columns
for seg in 1:n_seg
    df[!, Symbol("T_seg_$(seg)")] = [fs.tether_tension[seg] for fs in force_states]
end

# Add per-ring compression columns
for ring in 1:p.n_rings
    df[!, Symbol("C_ring_$(ring)")] = [fs.ring_compression[ring] for fs in force_states]
end
```

### Step 3: Add `include("force_analysis.jl")` to `main.jl`

`main.jl` already includes `parameters.jl`, `wind_profile.jl`, `dynamics.jl`,
`visualization.jl`. Add `include("force_analysis.jl")` and `include("geometry.jl")`
after `include("dynamics.jl")`.

### Step 4: Test headless run

```bash
julia --project=. src/main.jl
```

Expected: prints solve stats, opens window, writes CSV. Check that the CSV has the
new columns:

```bash
julia --project=. -e 'using CSV, DataFrames; df = CSV.read("output/simulation_10kw.csv", DataFrame); println(names(df))'
```

Expected: column names include `T_seg_1`, `C_ring_1`, `T_max_run`, `tau_aero`, etc.

### Step 5: Commit

```bash
git add src/main.jl
git commit -m "feat: extend CSV export with per-element force columns and FoS calibration data"
```

---

## Task 5: Update memory and verify final state

### Step 1: Run full test suite one final time

```bash
julia --project=. test/runtests.jl
```

Expected: all tests pass. Count should be > 168 (the original baseline).

### Step 2: Quick smoke test of 50 kW dashboard

```bash
julia --project=. scripts/run_dashboard_50kw.jl
```

Verify the scaled geometry (R=11.18m, longer shaft) renders correctly with blades and
force coloring at 50 kW parameters.

### Step 3: Commit memory update

Update `memory/MEMORY.md` to reflect the new module structure and function signatures.

---

## Completion Checklist

- [ ] `src/geometry.jl` created — inclined shaft, blades, ground plane
- [ ] `src/force_analysis.jl` created — `ForceState`, `element_forces`, `run_force_scan`
- [ ] `src/visualization.jl` refactored — 3 internal builders, force coloring, FoS colourbar, dynamic torque panel
- [ ] `src/TRPTKiteTurbineSimulator.jl` updated — includes both new modules
- [ ] `test/test_geometry.jl` updated — tests inclined shaft geometry correctly
- [ ] `test/test_force_analysis.jl` created — 4 testsets
- [ ] `test/runtests.jl` updated — includes force analysis tests
- [ ] `src/main.jl` updated — extended CSV export
- [ ] All dashboard scripts updated — include geometry.jl and force_analysis.jl
- [ ] Full test suite passes
- [ ] Manual smoke test confirms window opens with new layout
