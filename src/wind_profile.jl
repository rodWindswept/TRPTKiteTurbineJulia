# src/wind_profile.jl
# Wind speed as a function of altitude using the Hellmann power law.

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
