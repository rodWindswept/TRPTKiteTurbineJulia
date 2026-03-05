# src/visualization.jl
# GLMakie 3D visualization of the TRPT tensegrity structure.
# Strategy: pre-compute ODE trajectory, then animate at ≥30 FPS.

using GLMakie
include("geometry.jl")

# ── Scene construction ────────────────────────────────────────────────────────

"""
    build_trpt_scene(p, traj)

Build the GLMakie Figure and Observables for TRPT animation.

`traj` is a NamedTuple with fields:
    .t          : Vector{Float64} — time points (s)
    .alpha_tot  : Vector{Float64} — total twist (rad) at each time
    .omega      : Vector{Float64} — rotor speed (rad/s)
    .power_kw   : Vector{Float64} — ground power output (kW)
    .v_hub      : Float64 or Vector{Float64} — hub wind speed (m/s).
                  Scalar → treated as constant throughout the run.
                  Vector → one value per frame; displayed live in the HUD.

Returns `(fig, time_obs)` where `time_obs::Observable{Int}` is the frame index.
"""
function build_trpt_scene(p::SystemParams, traj)
    n_frames  = length(traj.t)
    v_hub_vec = isa(traj.v_hub, AbstractVector) ? traj.v_hub : fill(traj.v_hub, n_frames)

    fig = Figure(size=(1400, 800))

    # ── Left: 3D scene ────────────────────────────────────────────────────────
    ax3d = Axis3(fig[1, 1],
                 title  = "TRPT Kite Turbine — Live Geometry",
                 xlabel = "X (m)", ylabel = "Y (m)", zlabel = "Altitude (m)",
                 aspect = :data)

    # ── Right: telemetry HUD ──────────────────────────────────────────────────
    hud = GridLayout(fig[1, 2])
    Label(hud[1, 1], "Telemetry"; fontsize=16, font=:bold, halign=:left)

    time_label   = Label(hud[2, 1], "t = 0.00 s";        halign=:left)
    power_label  = Label(hud[3, 1], "P = 0.00 kW";       halign=:left)
    omega_label  = Label(hud[4, 1], "ω = 0.00 rad/s";    halign=:left)
    twist_label  = Label(hud[5, 1], "α_seg = 0.00°";     halign=:left)
    margin_label = Label(hud[6, 1], "Collapse margin: 100%"; halign=:left)
    wind_label   = Label(hud[7, 1], "V_hub = $(round(v_hub_vec[1], digits=2)) m/s"; halign=:left)
    Label(hud[8, 1], ""; halign=:left)   # spacer

    # Time slider
    Label(hud[9, 1], "Time"; halign=:left)
    time_slider = Slider(hud[10, 1], range=1:n_frames, startvalue=1)

    # Play / Pause button
    play_button = Button(hud[11, 1], label="▶ Play")

    colsize!(fig.layout, 2, Fixed(280))

    # ── Observable frame index ────────────────────────────────────────────────
    time_obs = time_slider.value

    # ── Reactive geometry: recomputes on every frame change ──────────────────
    nodes_obs = @lift begin
        compute_trpt_geometry(p, traj.alpha_tot[$time_obs])
    end

    # ── Draw tether lines (blue) ──────────────────────────────────────────────
    for j in 1:p.n_lines
        xs = @lift $nodes_obs[:, j, 1]
        ys = @lift $nodes_obs[:, j, 2]
        zs = @lift $nodes_obs[:, j, 3]
        lines!(ax3d, xs, ys, zs, color=:royalblue, linewidth=1.5)
    end

    # ── Draw polygon rings (black) ────────────────────────────────────────────
    n_levels = p.n_rings + 2
    for i in 2:(n_levels - 1)
        ring_x = @lift [$nodes_obs[i, :, 1]; $nodes_obs[i, 1, 1]]
        ring_y = @lift [$nodes_obs[i, :, 2]; $nodes_obs[i, 1, 2]]
        ring_z = @lift [$nodes_obs[i, :, 3]; $nodes_obs[i, 1, 3]]
        lines!(ax3d, ring_x, ring_y, ring_z, color=:black, linewidth=1.0)
    end

    # ── Draw rotor ring (red) ─────────────────────────────────────────────────
    rotor_x = @lift [$nodes_obs[end, :, 1]; $nodes_obs[end, 1, 1]]
    rotor_y = @lift [$nodes_obs[end, :, 2]; $nodes_obs[end, 1, 2]]
    rotor_z = @lift [$nodes_obs[end, :, 3]; $nodes_obs[end, 1, 3]]
    lines!(ax3d, rotor_x, rotor_y, rotor_z, color=:firebrick, linewidth=3.5)

    # ── Ground anchor ─────────────────────────────────────────────────────────
    scatter!(ax3d, [0.0], [0.0], [0.0], color=:green, markersize=20)

    # ── HUD update callback ───────────────────────────────────────────────────
    n_seg = p.n_rings + 1
    on(time_obs) do frame_idx
        t       = traj.t[frame_idx]
        alpha   = traj.alpha_tot[frame_idx]
        omega   = traj.omega[frame_idx]
        power   = traj.power_kw[frame_idx]
        alpha_s = rad2deg(alpha / n_seg)
        margin  = max(0.0, (1.0 - abs(alpha / n_seg) / π) * 100.0)

        time_label.text[]   = "t = $(round(t,     digits=2)) s"
        power_label.text[]  = "P = $(round(power,  digits=2)) kW"
        omega_label.text[]  = "ω = $(round(omega,  digits=3)) rad/s"
        twist_label.text[]  = "α_seg = $(round(alpha_s, digits=1))°"
        margin_label.text[] = "Collapse margin: $(round(margin, digits=1))%"
        wind_label.text[]   = "V_hub = $(round(v_hub_vec[frame_idx], digits=2)) m/s"
    end

    # ── Play button toggles animation ─────────────────────────────────────────
    is_playing = Observable(false)
    on(play_button.clicks) do _
        is_playing[] = !is_playing[]
        play_button.label[] = is_playing[] ? "⏸ Pause" : "▶ Play"
    end

    # Animation loop at ~30 FPS; stops at last frame
    @async while true
        if is_playing[]
            next_frame = min(time_obs[] + 1, n_frames)
            set_close_to!(time_slider, next_frame)
            if next_frame == n_frames
                is_playing[] = false
                play_button.label[] = "▶ Play"
            end
        end
        sleep(1 / 30)
    end

    return fig, time_obs
end
