# src/visualization.jl
# GLMakie 3D visualization of the TRPT tensegrity structure.
# Depends on: geometry.jl, force_analysis.jl, parameters.jl, dynamics.jl
# Strategy: pre-compute ODE trajectory + force scan, then animate at 30 FPS.

using GLMakie
using LinearAlgebra
using DifferentialEquations
using Printf

# ── Colour helper ──────────────────────────────────────────────────────────────

"""Map scalar `v` in [0, v_max] to blue→red RGBf colour."""
function _force_color(v::Float64, v_max::Float64)
    t = v_max > 0.0 ? clamp(v / v_max, 0.0, 1.0) : 0.0
    return RGBf(t, 0.0, 1.0 - t)
end

# ── 3D axes builder ────────────────────────────────────────────────────────────

function _build_3d_axes!(fig, position, p, p_obs, traj_obs, shaft_dir_obs,
                          T_max_obs, C_max_obs, time_obs)
    ax3d = Axis3(fig[position...],
                 title  = "TRPT Kite Turbine — Live Geometry",
                 xlabel = "X downwind (m)", ylabel = "Y crosswind (m)",
                 zlabel = "Altitude (m)",
                 aspect = :data)

    # Ground plane (static)
    for (xs, ys, zs) in world_ground_plane()
        lines!(ax3d, xs, ys, zs; color=(:grey, 0.3), linewidth=0.5)
    end

    # Wind direction arrow — updates with azimuth slider
    wind_xs = @lift begin
        sd = $shaft_dir_obs
        [0.0, sd[1] * 15.0]
    end
    wind_ys = @lift begin
        sd = $shaft_dir_obs
        [0.0, sd[2] * 15.0]
    end
    wind_zs = Observable([0.5, 0.5])
    lines!(ax3d, wind_xs, wind_ys, wind_zs; color=:darkorange, linewidth=3.0)

    # TRPT node positions — react to time AND shaft_dir
    nodes_obs = @lift begin
        compute_trpt_geometry(p, ($traj_obs).alpha_tot[$time_obs], $shaft_dir_obs)
    end

    # Blade positions
    blades_obs = @lift begin
        compute_blade_geometry(p, ($traj_obs).alpha_tot[$time_obs], $shaft_dir_obs)
    end

    n_seg    = p.n_rings + 1
    n_levels = p.n_rings + 2

    # Tether lines — force-coloured per segment
    for j in 1:p.n_lines
        xs = @lift $nodes_obs[:, j, 1]
        ys = @lift $nodes_obs[:, j, 2]
        zs = @lift $nodes_obs[:, j, 3]

        lc = @lift begin
            traj  = $traj_obs
            u     = [traj.alpha_tot[$time_obs], traj.omega[$time_obs]]
            vhv   = traj.v_hub
            v_hub = vhv isa AbstractVector ? vhv[$time_obs] : vhv
            fs    = element_forces($p_obs, u, v_hub)
            colors = Vector{RGBf}(undef, n_levels)
            colors[1] = _force_color(fs.tether_tension[1], TETHER_SWL)
            for ii in 2:n_levels
                seg_idx = min(ii - 1, n_seg)
                colors[ii] = _force_color(fs.tether_tension[seg_idx], TETHER_SWL)
            end
            colors
        end
        lines!(ax3d, xs, ys, zs; color=lc, linewidth=1.5)
    end

    # Ring polygons — force-coloured by compression
    for i in 2:(n_levels - 1)
        ring_x = @lift [$nodes_obs[i, :, 1]; $nodes_obs[i, 1, 1]]
        ring_y = @lift [$nodes_obs[i, :, 2]; $nodes_obs[i, 1, 2]]
        ring_z = @lift [$nodes_obs[i, :, 3]; $nodes_obs[i, 1, 3]]

        rc = @lift begin
            traj  = $traj_obs
            u     = [traj.alpha_tot[$time_obs], traj.omega[$time_obs]]
            vhv   = traj.v_hub
            v_hub = vhv isa AbstractVector ? vhv[$time_obs] : vhv
            fs    = element_forces($p_obs, u, v_hub)
            _force_color(fs.ring_compression[i - 1], RING_SWL)
        end
        lines!(ax3d, ring_x, ring_y, ring_z; color=rc, linewidth=1.0)
    end

    # Rotor ring (firebrick landmark)
    rotor_x = @lift [$nodes_obs[end, :, 1]; $nodes_obs[end, 1, 1]]
    rotor_y = @lift [$nodes_obs[end, :, 2]; $nodes_obs[end, 1, 2]]
    rotor_z = @lift [$nodes_obs[end, :, 3]; $nodes_obs[end, 1, 3]]
    lines!(ax3d, rotor_x, rotor_y, rotor_z; color=:firebrick, linewidth=3.5)

    # Blades
    for b in 1:p.n_blades
        bx = @lift $blades_obs[b, [1,2,3,4,1], 1]
        by = @lift $blades_obs[b, [1,2,3,4,1], 2]
        bz = @lift $blades_obs[b, [1,2,3,4,1], 3]
        lines!(ax3d, bx, by, bz; color=:steelblue, linewidth=2.0)
    end

    # Ground anchor
    scatter!(ax3d, [0.0], [0.0], [0.0]; color=:green, markersize=20)

    # Lifter line system — physically correct topology:
    #   Each top-ring node → common bearing point (on shaft axis, ~0.5 segments above rotor hub)
    #   Bearing point → lift point (~1 m further along shaft axis)
    #   At lift point: lifter kite force and anchor force meet
    l_seg          = p.tether_length / (p.n_rings + 1)
    bearing_offset = 1.5 * l_seg   # distance above rotor hub to bearing (m) — 3× segment spacing
    lift_offset    = 1.0            # distance above bearing to lift point (m)

    bearing_obs = @lift begin
        sd = $shaft_dir_obs
        rotor_hub = p.tether_length .* sd
        rotor_hub .+ bearing_offset .* sd
    end

    lift_point_obs = @lift begin
        bp = $bearing_obs
        sd = $shaft_dir_obs
        bp .+ lift_offset .* sd
    end

    # Lines from each top-ring node to the common bearing
    for j in 1:p.n_lines
        node_to_bearing_obs = @lift begin
            nd = $nodes_obs
            bp = $bearing_obs
            node = nd[end, j, :]
            ([node[1], bp[1]], [node[2], bp[2]], [node[3], bp[3]])
        end
        lines!(ax3d, @lift($node_to_bearing_obs[1]),
                     @lift($node_to_bearing_obs[2]),
                     @lift($node_to_bearing_obs[3]);
               color=:gold, linewidth=1.5)
    end

    # Single lifting line from bearing to lift point
    bearing_to_lift_obs = @lift begin
        bp = $bearing_obs
        lp = $lift_point_obs
        ([bp[1], lp[1]], [bp[2], lp[2]], [bp[3], lp[3]])
    end
    lines!(ax3d, @lift($bearing_to_lift_obs[1]),
                 @lift($bearing_to_lift_obs[2]),
                 @lift($bearing_to_lift_obs[3]);
           color=:gold, linewidth=3.0)

    # Markers at bearing and lift point
    scatter!(ax3d, @lift([$bearing_obs[1]]),
                   @lift([$bearing_obs[2]]),
                   @lift([$bearing_obs[3]]);
             color=:gold, markersize=10)
    scatter!(ax3d, @lift([$lift_point_obs[1]]),
                   @lift([$lift_point_obs[2]]),
                   @lift([$lift_point_obs[3]]);
             color=:white, markersize=8)

    return ax3d
end

# ── HUD builder ────────────────────────────────────────────────────────────────

function _build_hud!(layout, p, p_obs, traj_obs, time_obs,
                     T_max_obs, C_max_obs,
                     P_max_obs, omega_max_obs, v_max_obs, alpha_s_max_obs)
    # Fixed column width prevents label jitter as numbers change width
    colsize!(layout, 1, Fixed(320))

    lbl(row, txt; kw...) = Label(layout[row, 1], txt;
                                  halign=:left, tellwidth=false, justification=:left, kw...)

    lbl(1, "Live Telemetry"; fontsize=16, font=:bold)
    time_lbl   = lbl(2, "Time  t =     0.00 s")
    power_lbl  = lbl(3, "Output power  P =   0.00 kW")
    omega_lbl  = lbl(4, "Rotor speed  ω =   0.000 rad/s  (   0.0 rpm)")
    twist_lbl  = lbl(5, "Shaft twist / section  α =   0.0°")
    margin_lbl = lbl(6, "Collapse margin:  100.0%")
    wind_lbl   = lbl(7, "Wind at hub  V =   0.00 m/s")
    lbl(8, "")

    # ── Structural load panel with colour scale guide ──────────────────────
    lbl(9, "Structural Loads (this frame)"; fontsize=14, font=:bold)

    force_cmap = cgrad([RGBf(0.0, 0.0, 1.0), RGBf(1.0, 0.0, 0.0)])

    lbl(10, "Tether tension  (SWL = $(Int(TETHER_SWL)) N)";
        fontsize=12, font=:bold, color=:steelblue)
    t_minmax_lbl = lbl(11, "min     0 N  ·  max     0 N  ·  FoS  —")
    Colorbar(layout[12, 1]; colormap=force_cmap, limits=(0.0, Float64(TETHER_SWL)),
             vertical=false, height=16, tellheight=true, tellwidth=false,
             label="0 N (blue)  →  $(Int(TETHER_SWL)) N SWL (red)",
             labelsize=9, ticksize=4, ticklabelsize=8)

    lbl(13, "Ring compression  (SWL = $(Int(RING_SWL)) N)";
        fontsize=12, font=:bold, color=:firebrick)
    c_minmax_lbl = lbl(14, "min     0 N  ·  max     0 N  ·  FoS  —")
    Colorbar(layout[15, 1]; colormap=force_cmap, limits=(0.0, Float64(RING_SWL)),
             vertical=false, height=16, tellheight=true, tellwidth=false,
             label="0 N (blue)  →  $(Int(RING_SWL)) N SWL (red)",
             labelsize=9, ticksize=4, ticklabelsize=8)

    lbl(16, "")   # spacer

    # ── Run-wide peak summary ──────────────────────────────────────────────
    lbl(17, "Run peaks (all frames)"; fontsize=12, font=:bold)
    t_peak_lbl     = lbl(18, "T_peak       0 N  ·  FoS  —")
    c_peak_lbl     = lbl(19, "C_peak       0 N  ·  FoS  —")
    p_peak_lbl     = lbl(20, "P_peak    0.00 kW")
    omega_peak_lbl = lbl(21, "ω_peak   0.000 rad/s  (   0.0 rpm)")
    v_peak_lbl     = lbl(22, "V_peak    0.00 m/s")
    twist_peak_lbl = lbl(23, "α_peak     0.0° / section")

    n_seg = p.n_rings + 1
    fos_str(v) = (isinf(v) || isnan(v)) ? "  ∞" : @sprintf("%5.1f", v)

    on(time_obs) do fi
        traj    = traj_obs[]
        t       = traj.t[fi]
        alpha   = traj.alpha_tot[fi]
        omega   = traj.omega[fi]
        power   = traj.power_kw[fi]
        vhv     = traj.v_hub
        v_hub   = vhv isa AbstractVector ? vhv[fi] : vhv
        alpha_s = rad2deg(alpha / n_seg)
        margin  = max(0.0, (1.0 - abs(alpha / n_seg) / π) * 100.0)
        rpm     = omega * 60.0 / (2π)

        time_lbl.text[]   = @sprintf("Time  t = %8.2f s",                    t)
        power_lbl.text[]  = @sprintf("Output power  P = %6.2f kW",           power)
        omega_lbl.text[]  = @sprintf("Rotor speed  ω = %7.3f rad/s  (%6.1f rpm)", omega, rpm)
        twist_lbl.text[]  = @sprintf("Shaft twist / section  α = %5.1f°",    alpha_s)
        margin_lbl.text[] = @sprintf("Collapse margin: %5.1f%%",             margin)
        wind_lbl.text[]   = @sprintf("Wind at hub  V = %5.2f m/s",           v_hub)

        u_frame = [alpha, omega]
        fs      = element_forces(p_obs[], u_frame, v_hub)
        T_min_f = minimum(fs.tether_tension)
        T_max_f = maximum(fs.tether_tension)
        C_min_f = minimum(fs.ring_compression)
        C_max_f = maximum(fs.ring_compression)
        T_fos_f = T_max_f > 0.0 ? TETHER_SWL / T_max_f : Inf
        C_fos_f = C_max_f > 0.0 ? RING_SWL   / C_max_f : Inf

        t_minmax_lbl.text[] = @sprintf("min %5.0f N  ·  max %5.0f N  ·  FoS %s",
                                        T_min_f, T_max_f, fos_str(T_fos_f))
        c_minmax_lbl.text[] = @sprintf("min %5.0f N  ·  max %5.0f N  ·  FoS %s",
                                        C_min_f, C_max_f, fos_str(C_fos_f))

        T_peak     = T_max_obs[]
        C_peak     = C_max_obs[]
        P_peak     = P_max_obs[]
        omega_peak = omega_max_obs[]
        v_peak     = v_max_obs[]
        alpha_s_pk = alpha_s_max_obs[]
        rpm_peak   = omega_peak * 60.0 / (2π)

        t_peak_lbl.text[]     = @sprintf("T_peak  %5.0f N  ·  FoS %s",
                                          T_peak, fos_str(T_peak > 0 ? TETHER_SWL / T_peak : Inf))
        c_peak_lbl.text[]     = @sprintf("C_peak  %5.0f N  ·  FoS %s",
                                          C_peak, fos_str(C_peak > 0 ? RING_SWL   / C_peak : Inf))
        p_peak_lbl.text[]     = @sprintf("P_peak  %6.2f kW", P_peak)
        omega_peak_lbl.text[] = @sprintf("ω_peak  %7.3f rad/s  (%6.1f rpm)", omega_peak, rpm_peak)
        v_peak_lbl.text[]     = @sprintf("V_peak  %5.2f m/s", v_peak)
        twist_peak_lbl.text[] = @sprintf("α_peak  %5.1f° / section", rad2deg(alpha_s_pk))
    end
end

# ── Controls builder ───────────────────────────────────────────────────────────

function _build_controls!(layout, p, p_obs, traj_obs,
                          T_max_obs, C_max_obs,
                          P_max_obs, omega_max_obs, v_max_obs, alpha_s_max_obs,
                          n_frames)
    Label(layout[1, 1], "Controls"; fontsize=14, font=:bold, halign=:left)

    elev_val_lbl = Label(layout[2, 1],
                         "Shaft tilt above horizon  β = $(round(rad2deg(p.elevation_angle), digits=0))°";
                         halign=:left)
    elev_slider = Slider(layout[3, 1]; range=23.0:1.0:75.0,
                         startvalue=clamp(rad2deg(p.elevation_angle), 23.0, 75.0))
    on(elev_slider.value) do v
        elev_val_lbl.text[] = "Shaft tilt above horizon  β = $(round(v, digits=0))°"
    end

    azimuth_val_lbl = Label(layout[4, 1],
                             "Wind direction  φ = 0°"; halign=:left)
    azimuth_slider = Slider(layout[5, 1]; range=0.0:1.0:360.0, startvalue=0.0)
    Label(layout[6, 1], "Visual only — yaws geometry; physics always head-to-wind";
          halign=:left, fontsize=9, color=:grey60)
    on(azimuth_slider.value) do v
        azimuth_val_lbl.text[] = "Wind direction  φ = $(round(v, digits=0))°"
    end

    shaft_dir_obs = Observable([cos(p.elevation_angle), 0.0, sin(p.elevation_angle)])

    update_shaft! = _ -> begin
        β  = deg2rad(elev_slider.value[])
        φw = deg2rad(azimuth_slider.value[])
        shaft_dir_obs[] = [cos(φw)*cos(β), sin(φw)*cos(β), sin(β)]
    end
    on(update_shaft!, elev_slider.value)
    on(update_shaft!, azimuth_slider.value)

    Label(layout[7, 1], "Time"; halign=:left)
    time_slider = Slider(layout[8, 1]; range=1:n_frames, startvalue=1)

    play_row    = GridLayout(layout[9, 1])
    play_btn    = Button(play_row[1, 1]; label="▶ Play")
    Label(play_row[1, 2], "Solver:"; halign=:right)
    solver_menu = Menu(play_row[1, 3]; options=["Tsit5", "RK4", "Euler"], default="Tsit5")
    Label(play_row[2, 1:3],
          "Tsit5 — adaptive  ·  RK4 — fixed-step  ·  Euler — fast, may drift";
          fontsize=9, halign=:left, color=:grey60)

    is_playing = Observable(false)
    on(play_btn.clicks) do _
        is_playing[] = !is_playing[]
        play_btn.label[] = is_playing[] ? "⏸ Pause" : "▶ Play"
    end
    @async while true
        if is_playing[]
            nf = min(time_slider.value[] + 1, n_frames)
            set_close_to!(time_slider, nf)
            if nf == n_frames
                is_playing[] = false
                play_btn.label[] = "▶ Play"
            end
        end
        sleep(1 / 30)
    end

    # Dynamic simulation re-run panel
    Label(layout[10, 1], "Simulation Re-run"; fontsize=13, font=:bold, halign=:left)
    enable_toggle = Toggle(layout[11, 1])
    Label(layout[11, 2], "Unlock re-run"; halign=:left)
    Label(layout[12, 1], "Replaces trajectory · uses β, V_ref, k_mppt below";
          halign=:left, fontsize=9, color=:grey60)

    Label(layout[13, 1], "Wind speed  V_ref (m/s)"; halign=:left)
    v_wind_slider = Slider(layout[14, 1]; range=4.0:0.5:20.0, startvalue=p.v_wind_ref)
    v_wind_lbl = Label(layout[15, 1],
                       @sprintf("%.1f m/s", p.v_wind_ref); halign=:left)
    on(v_wind_slider.value) do v
        v_wind_lbl.text[] = @sprintf("%.1f m/s", v)
    end

    Label(layout[16, 1], "MPPT gain  k_mppt (N·m·s²/rad²)"; halign=:left)
    k_mppt_range = vcat([0.0], exp10.(range(-1.0, log10(5000.0), length=199)))
    k_mppt_slider = Slider(layout[17, 1]; range=k_mppt_range, startvalue=p.k_mppt)
    k_mppt_lbl = Label(layout[18, 1],
                       @sprintf("%.2f N·m·s²/rad²", p.k_mppt); halign=:left)
    on(k_mppt_slider.value) do v
        k_mppt_lbl.text[] = v == 0.0 ? "0  (freewheel)" :
                                        @sprintf("%.2f N·m·s²/rad²", v)
    end

    rerun_btn = Button(layout[19, 1]; label="Re-run ODE (120 s)")
    on(rerun_btn.clicks) do _
        enable_toggle.active[] || return

        k_new   = k_mppt_slider.value[]
        v_new   = v_wind_slider.value[]
        β_new   = deg2rad(elev_slider.value[])   # shaft tilt from elevation slider
        p_new = SystemParams(p.rho, v_new, p.h_ref, β_new,
                             p.lifter_elevation,
                             p.rotor_radius, p.tether_length, p.trpt_hub_radius,
                             p.trpt_rL_ratio, p.n_lines, p.tether_diameter,
                             p.e_modulus, p.n_rings, p.m_ring, p.n_blades,
                             p.m_blade, p.cp, p.i_pto, k_new)

        solver = if solver_menu.selection[] == "Tsit5"
            Tsit5()
        elseif solver_menu.selection[] == "RK4"
            RK4()
        else
            Euler()
        end

        # Euler is fixed-step — reltol/abstol do not apply; adaptive solvers use both
        kwargs = solver_menu.selection[] == "Euler" ?
                 (saveat=1/30, dt=0.01) :
                 (reltol=1e-6, abstol=1e-6, saveat=1/30)

        sol_new = solve(ODEProblem(trpt_ode!, [0.0, 1.0], (0.0, 120.0), p_new),
                        solver; kwargs...)

        h_hub_new   = hub_altitude(p_new.tether_length, p_new.elevation_angle)
        v_hub_new   = wind_at_altitude(p_new.v_wind_ref, p_new.h_ref, h_hub_new)
        n_frames_new = length(sol_new.t)
        # v_hub must be Vector{Float64} to match the type of traj_obs (set at scene build)
        v_hub_vec = fill(v_hub_new, n_frames_new)
        new_traj = (
            t         = sol_new.t,
            alpha_tot = sol_new[1, :],
            omega     = sol_new[2, :],
            power_kw  = [instantaneous_power(p_new, [sol_new[1,i], sol_new[2,i]]) / 1000.0
                         for i in eachindex(sol_new.t)],
            v_hub     = v_hub_vec,
        )

        u_frames = [[new_traj.alpha_tot[i], new_traj.omega[i]] for i in eachindex(new_traj.t)]
        T_new, C_new = run_force_scan(p_new, u_frames, v_hub_vec)

        n_seg_new = p_new.n_rings + 1
        p_obs[]           = p_new
        traj_obs[]        = new_traj
        T_max_obs[]       = T_new
        C_max_obs[]       = C_new
        P_max_obs[]       = maximum(new_traj.power_kw)
        omega_max_obs[]   = maximum(abs.(new_traj.omega))
        v_max_obs[]       = v_hub_new
        alpha_s_max_obs[] = maximum(abs.(new_traj.alpha_tot)) / n_seg_new
        set_close_to!(time_slider, 1)
        # Force HUD refresh — set_close_to! is a no-op when slider is already at 1,
        # so the on(time_obs) callback would not fire and peak labels would stay stale.
        notify(time_slider.value)
    end

    return time_slider, shaft_dir_obs
end

# ── Main scene entry point ─────────────────────────────────────────────────────

"""
    build_trpt_scene(p, traj)

Build the GLMakie Figure for interactive TRPT visualization.

`traj` fields: `.t`, `.alpha_tot`, `.omega`, `.power_kw`, `.v_hub` (scalar or vector).
Returns `(fig, time_obs)`.
"""
function build_trpt_scene(p::SystemParams, traj)
    n_frames  = length(traj.t)
    vhv       = isa(traj.v_hub, AbstractVector) ? traj.v_hub : fill(traj.v_hub, n_frames)

    u_frames  = [[traj.alpha_tot[i], traj.omega[i]] for i in 1:n_frames]
    T_max_run, C_max_run = run_force_scan(p, u_frames, vhv)

    n_seg     = p.n_rings + 1
    traj_norm = (t=traj.t, alpha_tot=traj.alpha_tot, omega=traj.omega,
                 power_kw=traj.power_kw, v_hub=vhv)
    p_obs           = Observable(p)
    traj_obs        = Observable(traj_norm)
    T_max_obs       = Observable(T_max_run)
    C_max_obs       = Observable(C_max_run)
    P_max_obs       = Observable(maximum(traj_norm.power_kw))
    omega_max_obs   = Observable(maximum(abs.(traj_norm.omega)))
    v_max_obs       = Observable(maximum(vhv))
    alpha_s_max_obs = Observable(maximum(abs.(traj_norm.alpha_tot)) / n_seg)

    fig = Figure(size=(1600, 900))

    right = GridLayout(fig[1, 2])
    colsize!(fig.layout, 2, Fixed(660))
    hud_l = GridLayout(right[1, 1])
    ctl_l = GridLayout(right[1, 2])
    colsize!(right, 1, Fixed(330))
    colsize!(right, 2, Fixed(330))

    time_slider, shaft_dir_obs = _build_controls!(ctl_l, p, p_obs, traj_obs,
                                                   T_max_obs, C_max_obs,
                                                   P_max_obs, omega_max_obs,
                                                   v_max_obs, alpha_s_max_obs,
                                                   n_frames)
    time_obs = time_slider.value

    _build_3d_axes!(fig, (1, 1), p, p_obs, traj_obs, shaft_dir_obs,
                    T_max_obs, C_max_obs, time_obs)
    _build_hud!(hud_l, p, p_obs, traj_obs, time_obs,
                T_max_obs, C_max_obs,
                P_max_obs, omega_max_obs, v_max_obs, alpha_s_max_obs)

    # Fire initial HUD update at frame 1 (callbacks registered after slider creation)
    notify(time_obs)

    return fig, time_obs
end
