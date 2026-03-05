# src/visualization.jl
# GLMakie 3D visualization of the TRPT tensegrity structure.
# Depends on: geometry.jl, force_analysis.jl, parameters.jl, dynamics.jl
# Strategy: pre-compute ODE trajectory + force scan, then animate at 30 FPS.

using GLMakie
using LinearAlgebra
using DifferentialEquations

# ── Colour helper ──────────────────────────────────────────────────────────────

"""Map scalar `v` in [0, v_max] to blue→red RGBf colour."""
function _force_color(v::Float64, v_max::Float64)
    t = v_max > 0.0 ? clamp(v / v_max, 0.0, 1.0) : 0.0
    return RGBf(t, 0.0, 1.0 - t)
end

# ── 3D axes builder ────────────────────────────────────────────────────────────

function _build_3d_axes!(fig, position, p, traj_obs, shaft_dir_obs,
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
            fs    = element_forces(p, u, v_hub)
            T_max = $T_max_obs
            colors = Vector{RGBf}(undef, n_levels)
            colors[1] = _force_color(fs.tether_tension[1], T_max)
            for ii in 2:n_levels
                seg_idx = min(ii - 1, n_seg)
                colors[ii] = _force_color(fs.tether_tension[seg_idx], T_max)
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
            fs    = element_forces(p, u, v_hub)
            C_max = $C_max_obs
            _force_color(fs.ring_compression[i - 1], C_max)
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

    # Lifter kite force arrows — one per top-ring node, showing lifter line direction
    # Arrow direction: same horizontal azimuth as shaft, at lifter_elevation above horizontal
    m_airborne_vis = p.n_blades * p.m_blade + p.n_rings * p.m_ring
    T_lift_vis     = m_airborne_vis * 9.81 / sin(p.lifter_elevation)
    arrow_len      = 5.0   # fixed visual length (m) — not force-scaled, just directional indicator

    for j in 1:p.n_lines
        ax_obs = @lift begin
            sd      = $shaft_dir_obs
            nd      = $nodes_obs
            # Horizontal direction of shaft (unit vector in ground plane)
            sd_horiz_mag = sqrt(sd[1]^2 + sd[2]^2)
            horiz = sd_horiz_mag > 1e-6 ? [sd[1]/sd_horiz_mag, sd[2]/sd_horiz_mag, 0.0] :
                                           [1.0, 0.0, 0.0]
            # Lifter line direction vector
            lift_dir = horiz .* cos(p.lifter_elevation) .+ [0.0, 0.0, sin(p.lifter_elevation)]
            node_pos = nd[end, j, :]
            tip_pos  = node_pos .+ arrow_len .* lift_dir
            ([node_pos[1], tip_pos[1]], [node_pos[2], tip_pos[2]], [node_pos[3], tip_pos[3]])
        end
        lines!(ax3d, @lift($ax_obs[1]), @lift($ax_obs[2]), @lift($ax_obs[3]);
               color=:gold, linewidth=2.0)
    end

    return ax3d
end

# ── HUD builder ────────────────────────────────────────────────────────────────

function _build_hud!(layout, p, traj_obs, time_obs, T_max_obs, C_max_obs)
    Label(layout[1, 1], "Live Telemetry"; fontsize=16, font=:bold, halign=:left)
    time_lbl   = Label(layout[2, 1], "Time  t = 0.00 s";                    halign=:left)
    power_lbl  = Label(layout[3, 1], "Output power  P = 0.00 kW";           halign=:left)
    omega_lbl  = Label(layout[4, 1], "Rotor speed  ω = 0.00 rad/s (0 rpm)"; halign=:left)
    twist_lbl  = Label(layout[5, 1], "Shaft twist / section  α = 0.00°";    halign=:left)
    margin_lbl = Label(layout[6, 1], "Collapse margin: 100%";               halign=:left)
    wind_lbl   = Label(layout[7, 1], "Wind at hub  V = 0.00 m/s";           halign=:left)
    Label(layout[8, 1], ""; halign=:left)

    Label(layout[9,  1], "Tether pull (axial tension, N)";    fontsize=12, font=:bold, halign=:left)
    t_bar_lbl = Label(layout[10, 1], "0 → — N  |  FoS: —"; halign=:left)
    Label(layout[11, 1], "Ring squeeze (hoop compression, N)"; fontsize=12, font=:bold, halign=:left)
    c_bar_lbl = Label(layout[12, 1], "0 → — N  |  FoS: —"; halign=:left)

    n_seg = p.n_rings + 1

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

        rpm = omega * 60.0 / (2π)
        time_lbl.text[]   = "Time  t = $(round(t,      digits=2)) s"
        power_lbl.text[]  = "Output power  P = $(round(power,  digits=2)) kW"
        omega_lbl.text[]  = "Rotor speed  ω = $(round(omega, digits=3)) rad/s  ($(round(rpm, digits=1)) rpm)"
        twist_lbl.text[]  = "Shaft twist / section  α = $(round(alpha_s, digits=1))°"
        margin_lbl.text[] = "Collapse margin: $(round(margin, digits=1))%"
        wind_lbl.text[]   = "Wind at hub  V = $(round(v_hub, digits=2)) m/s"

        T_max = T_max_obs[]
        C_max = C_max_obs[]
        T_fos = T_max > 0 ? round(TETHER_SWL / T_max, digits=1) : Inf
        C_fos = C_max > 0 ? round(RING_SWL   / C_max, digits=1) : Inf
        t_bar_lbl.text[] = "0 → $(round(T_max, digits=1)) N  |  SWL=$(TETHER_SWL)N  FoS=$(T_fos)"
        c_bar_lbl.text[] = "0 → $(round(C_max, digits=1)) N  |  SWL=$(RING_SWL)N  FoS=$(C_fos)"
    end
end

# ── Controls builder ───────────────────────────────────────────────────────────

function _build_controls!(layout, p, traj_obs, T_max_obs, C_max_obs, n_frames)
    Label(layout[1, 1], "Controls"; fontsize=14, font=:bold, halign=:left)

    Label(layout[2, 1], "Shaft tilt above horizon  β (°)"; halign=:left)
    elev_slider = Slider(layout[3, 1]; range=0.0:1.0:75.0,
                         startvalue=rad2deg(p.elevation_angle))

    Label(layout[4, 1], "Wind direction  φ (°)"; halign=:left)
    azimuth_slider = Slider(layout[5, 1]; range=0.0:1.0:360.0, startvalue=0.0)

    shaft_dir_obs = Observable([cos(p.elevation_angle), 0.0, sin(p.elevation_angle)])

    update_shaft! = _ -> begin
        β  = deg2rad(elev_slider.value[])
        φw = deg2rad(azimuth_slider.value[])
        shaft_dir_obs[] = [cos(φw)*cos(β), sin(φw)*cos(β), sin(β)]
    end
    on(update_shaft!, elev_slider.value)
    on(update_shaft!, azimuth_slider.value)

    Label(layout[6, 1], "Time"; halign=:left)
    time_slider = Slider(layout[7, 1]; range=1:n_frames, startvalue=1)

    play_row    = GridLayout(layout[8, 1])
    play_btn    = Button(play_row[1, 1]; label="▶ Play")
    Label(play_row[1, 2], "Solver:"; halign=:right)
    solver_menu = Menu(play_row[1, 3]; options=["Tsit5", "RK4", "Euler"], default="Tsit5")

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

    # Dynamic torque panel
    Label(layout[9, 1], "Power Take-Off Tuning"; fontsize=13, font=:bold, halign=:left)
    enable_toggle = Toggle(layout[10, 1])
    Label(layout[10, 2], "Enable"; halign=:left)

    Label(layout[11, 1], "Generator braking  c_pto (N·m·s/rad)"; halign=:left)
    c_pto_slider = Slider(layout[12, 1];
                          range=exp10.(range(log10(100.0), log10(50000.0), length=200)),
                          startvalue=p.c_pto)
    c_pto_lbl = Label(layout[13, 1], "$(round(p.c_pto, digits=0)) N·m·s/rad"; halign=:left)
    on(c_pto_slider.value) do v
        c_pto_lbl.text[] = "$(round(v, digits=0)) N·m·s/rad"
    end

    rerun_btn = Button(layout[14, 1]; label="Re-run ODE")
    on(rerun_btn.clicks) do _
        enable_toggle.active[] || return

        c_new = c_pto_slider.value[]
        p_new = SystemParams(p.rho, p.v_wind_ref, p.h_ref, p.elevation_angle,
                             p.lifter_elevation,
                             p.rotor_radius, p.tether_length, p.trpt_hub_radius,
                             p.trpt_rL_ratio, p.n_lines, p.tether_diameter,
                             p.e_modulus, p.n_rings, p.m_ring, p.n_blades,
                             p.m_blade, p.cp, p.i_pto, c_new)

        solver = if solver_menu.selection[] == "Tsit5"
            Tsit5()
        elseif solver_menu.selection[] == "RK4"
            RK4()
        else
            Euler()
        end

        kwargs = solver_menu.selection[] == "Euler" ?
                 (reltol=1e-6, abstol=1e-6, saveat=1/30, dt=0.01) :
                 (reltol=1e-6, abstol=1e-6, saveat=1/30)

        sol_new = solve(ODEProblem(trpt_ode!, [0.0, 1.0], (0.0, 120.0), p_new),
                        solver; kwargs...)

        old_traj = traj_obs[]
        new_traj = (
            t         = sol_new.t,
            alpha_tot = sol_new[1, :],
            omega     = sol_new[2, :],
            power_kw  = [instantaneous_power(p_new, [sol_new[1,i], sol_new[2,i]]) / 1000.0
                         for i in eachindex(sol_new.t)],
            v_hub     = old_traj.v_hub,
        )

        vhv = new_traj.v_hub
        v_hub_scalar = vhv isa AbstractVector ? vhv[1] : vhv
        u_frames  = [[new_traj.alpha_tot[i], new_traj.omega[i]] for i in eachindex(new_traj.t)]
        v_hubs    = fill(v_hub_scalar, length(new_traj.t))
        T_new, C_new = run_force_scan(p_new, u_frames, v_hubs)

        traj_obs[]  = new_traj
        T_max_obs[] = T_new
        C_max_obs[] = C_new
        set_close_to!(time_slider, 1)
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

    traj_norm = (t=traj.t, alpha_tot=traj.alpha_tot, omega=traj.omega,
                 power_kw=traj.power_kw, v_hub=vhv)
    traj_obs  = Observable(traj_norm)
    T_max_obs = Observable(T_max_run)
    C_max_obs = Observable(C_max_run)

    fig = Figure(size=(1600, 900))

    right = GridLayout(fig[1, 2])
    colsize!(fig.layout, 2, Fixed(340))
    hud_l = GridLayout(right[1, 1])
    ctl_l = GridLayout(right[2, 1])

    time_slider, shaft_dir_obs = _build_controls!(ctl_l, p, traj_obs,
                                                   T_max_obs, C_max_obs, n_frames)
    time_obs = time_slider.value

    _build_3d_axes!(fig, (1, 1), p, traj_obs, shaft_dir_obs,
                    T_max_obs, C_max_obs, time_obs)
    _build_hud!(hud_l, p, traj_obs, time_obs, T_max_obs, C_max_obs)

    return fig, time_obs
end
