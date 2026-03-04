# src/wind_profile.jl
# Wind speed as a function of altitude (Hellmann power law) and
# time-varying wind model constructors for simulation scenarios.

using Random

"""
    wind_at_altitude(v_ref, h_ref, h; hellmann_exponent = 1/7) -> Float64

Return the wind speed at altitude `h` using the Hellmann (power-law) wind profile:

    V(h) = V_ref × (h / h_ref)^α

where α is the Hellmann exponent (default 1/7, appropriate for open flat terrain).
Returns 0.0 for any non-positive altitude.

# Arguments
- `v_ref`             : Reference wind speed (m/s) measured at `h_ref`
- `h_ref`             : Reference altitude (m)
- `h`                 : Target altitude (m)
- `hellmann_exponent` : Terrain-dependent shear exponent α (dimensionless); default 1/7
"""
function wind_at_altitude(v_ref::Float64, h_ref::Float64, h::Float64;
                          hellmann_exponent::Float64 = 1.0 / 7.0)::Float64
    if h <= 0.0
        return 0.0
    end
    return v_ref * (h / h_ref)^hellmann_exponent
end

"""
    hub_altitude(tether_length, elevation_angle) -> Float64

Return the vertical altitude of the airborne rotor hub (m), computed from the
TRPT tether geometry as:

    h_hub = tether_length × sin(elevation_angle)

# Arguments
- `tether_length`   : Total TRPT tether length L₀ (m)
- `elevation_angle` : Shaft elevation angle β above horizontal (rad)
"""
function hub_altitude(tether_length::Float64, elevation_angle::Float64)::Float64
    return tether_length * sin(elevation_angle)
end

# ── Wind model constructors ────────────────────────────────────────────────────
# Each returns a closure  f(t::Float64)::Float64  giving v_ref (m/s at reference
# height p.h_ref). Pass the closure as `wind_fn` to trpt_ode_wind!.

"""
    steady_wind(v_ref) -> Function

Return a constant wind function: `v(t) = v_ref` for all t.
"""
steady_wind(v_ref::Float64) = (t::Float64) -> v_ref

"""
    wind_ramp(v_start, v_end, t_ramp_start, t_ramp_end) -> Function

Return a wind function that linearly ramps from `v_start` to `v_end` between
`t_ramp_start` and `t_ramp_end` (s). Constant outside that window.
"""
function wind_ramp(v_start::Float64, v_end::Float64,
                   t_ramp_start::Float64, t_ramp_end::Float64)
    function f(t::Float64)::Float64
        t <= t_ramp_start && return v_start
        t >= t_ramp_end   && return v_end
        return v_start + (v_end - v_start) * (t - t_ramp_start) / (t_ramp_end - t_ramp_start)
    end
    return f
end

"""
    gust_event(v_base, v_gust, t_start, t_end) -> Function

Return a wind function with a raised-cosine (Hann-window) gust from `t_start`
to `t_end`. Peak speed is `v_gust`; baseline speed is `v_base` outside the
gust window. Smooth (C¹ continuous) at the gust edges.
"""
function gust_event(v_base::Float64, v_gust::Float64,
                    t_start::Float64, t_end::Float64)
    function f(t::Float64)::Float64
        (t < t_start || t > t_end) && return v_base
        frac = (t - t_start) / (t_end - t_start)
        return v_base + (v_gust - v_base) * 0.5 * (1.0 - cos(2π * frac))
    end
    return f
end

"""
    turbulent_wind(v_mean, turbulence_intensity, t_max;
                   dt = 1/30, rng_seed = 42) -> Function

Return a wind function based on a first-order Markov (AR(1)) turbulence model.

- `turbulence_intensity`: σ/μ, e.g. 0.15 for 15 % TI (typical onshore Class A).
- Integral length scale: L = 340 m (IEC 61400-1 Class A at 30 m hub).
- Pre-computes a time series on `[0, t_max]` at step `dt`, then interpolates
  linearly for any query time.
- Wind speed is clamped to ≥ 0.5 m/s to prevent negative values.
"""
function turbulent_wind(v_mean::Float64, turbulence_intensity::Float64,
                        t_max::Float64;
                        dt::Float64   = 1.0 / 30.0,
                        rng_seed::Int = 42)
    σ   = turbulence_intensity * v_mean
    L   = 340.0               # IEC 61400-1 integral length scale (m)
    T_L = L / v_mean          # integral time scale (s)
    φ   = exp(-dt / T_L)      # AR(1) autocorrelation coefficient

    rng  = MersenneTwister(rng_seed)
    n    = round(Int, t_max / dt) + 2
    ts   = [(i - 1) * dt for i in 1:n]
    vs   = zeros(Float64, n)
    vs[1] = v_mean
    w     = 0.0
    for i in 2:n
        w     = φ * w + sqrt(1.0 - φ^2) * randn(rng)
        vs[i] = max(0.5, v_mean + σ * w)
    end

    function interp(t::Float64)::Float64
        i = searchsortedfirst(ts, t)
        i == 1             && return vs[1]
        i > length(ts)     && return vs[end]
        frac = (t - ts[i - 1]) / (ts[i] - ts[i - 1])
        return vs[i - 1] + frac * (vs[i] - vs[i - 1])
    end
    return interp
end
