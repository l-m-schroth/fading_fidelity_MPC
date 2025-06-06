import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from acados_template import latexify_plot
import os
from utils_shared import get_dir

def plot_drone_mpc_solution(
    mpc, 
    reference_xy=None,
    closed_loop_traj=None,
    open_loop_plan=None,
    u_traj=None,
    plot_title="Drone + Pendulum MPC",
    step_pose=10,
    number=None,        
    legend=True,        
    initial_state=None,  
    save=False,     
    latexify=True      
):
    """
    Updated plotting function for the new DroneMPC class logic.
    
    (Note: Only modifications in the configuration plot and addition of a fourth plot for r(t) 
     are made per the explicit instructions.)
    """

    if latexify:
        latexify_plot(fontsize=16)

    # ------------------------------------------------
    # 1) Decide data source: multi-phase open_loop_plan
    #    vs single-phase closed_loop_traj
    # ------------------------------------------------
    if open_loop_plan is not None:
        data_list = open_loop_plan
        multi_phase = True
    elif closed_loop_traj is not None:
        data_list = [row for row in closed_loop_traj]
        multi_phase = False
    else:
        print("No data to plot (both open_loop_plan and closed_loop_traj are None).")
        return

    if not data_list:
        print("Empty data list => nothing to plot.")
        return

    # ------------------------------------------------
    # Helpers to unify snapshots & controls
    # ------------------------------------------------
    def unify_snapshots_into_2d(snapshot_list):
        filtered = []
        for arr in snapshot_list:
            if arr is None:
                continue
            arr2 = np.asarray(arr)
            if arr2.size == 0:
                continue
            if arr2.ndim == 1:
                arr2 = arr2.reshape(1, -1)
            filtered.append(arr2)
        if len(filtered) == 0:
            return None
        return np.vstack(filtered)

    def unify_controls(u_list):
        filtered = []
        for item in u_list:
            if item is None:
                continue
            arr = np.asarray(item)
            if arr.size == 0:
                continue
            if arr.ndim == 1:
                arr = arr.reshape(1, -1)
            filtered.append(arr)
        if len(filtered)==0:
            return None
        return np.vstack(filtered)

    # We'll figure out if we have up to 3 phases for multi-phase:
    N = mpc.opts.N
    switch_stage = np.clip(mpc.opts.switch_stage, 0, N)

    # ------------------------------------------------
    # 2) Partition states for multi-phase vs single-phase
    # ------------------------------------------------
    if multi_phase:
        phaseA_list = data_list[:switch_stage+1]
        phaseB_list = data_list[switch_stage+1:]
        phaseA_data = unify_snapshots_into_2d(phaseA_list)
        phaseB_data = unify_snapshots_into_2d(phaseB_list)

        phases_data   = []
        phases_colors = []
        if phaseA_data is not None and phaseA_data.shape[0]>0:
            phases_data.append(phaseA_data)
            phases_colors.append("red")
        if phaseB_data is not None and phaseB_data.shape[0]>0:
            phases_data.append(phaseB_data)
            phases_colors.append("blue")

        phases_u = []
        phases_u_colors = []
        if u_traj is not None:
            switch_u = min(switch_stage, len(u_traj))
            uA_list = u_traj[:switch_u]
            uB_list = u_traj[switch_u:]
            phaseA_u = unify_controls(uA_list)
            phaseB_u = unify_controls(uB_list)

            if phaseA_data is not None and phaseA_data.shape[0]>0:
                phases_u.append(phaseA_u)
                phases_u_colors.append("red")
            if phaseB_data is not None and phaseB_data.shape[0]>0:
                phases_u.append(phaseB_u)
                phases_u_colors.append("blue")
        else:
            phases_u = []
            phases_u_colors = []

    else:
        single_data = unify_snapshots_into_2d(data_list)
        phases_data   = []
        phases_colors = []
        if single_data is not None and single_data.shape[0]>0:
            phases_data.append(single_data)
            phases_colors.append("red")

        phases_u = []
        phases_u_colors= []
        if u_traj is not None:
            single_u = unify_controls(u_traj)
            phases_u.append(single_u)
            phases_u_colors.append("red")
        else:
            phases_u = [None]
            phases_u_colors= ["red"]

    # ------------------------------------------------
    # FIGURE 1: (y,z) + time subplots
    # ------------------------------------------------
    fig1 = plt.figure(figsize=(10,10))
    gs = fig1.add_gridspec(4, 2, height_ratios=[1.5,1,1,1])
    ax_xy = fig1.add_subplot(gs[0,:])
    ax_xy.set_title(f"{plot_title} - YZ Trajectory")
    ax_xy.set_xlabel("y [m]")
    ax_xy.set_ylabel("z [m]")
    ax_xy.grid(True)
    ax_xy.axis("equal")

    if reference_xy is not None:
        ax_xy.plot(reference_xy[0], reference_xy[1], marker='*', color='gold', markersize=10, label="Target")

    ax_y   = fig1.add_subplot(gs[1,0]); ax_y.set_title("y(t)");   ax_y.grid(True)
    ax_z   = fig1.add_subplot(gs[1,1]); ax_z.set_title("z(t)");   ax_z.grid(True)
    ax_phi = fig1.add_subplot(gs[2,0]); ax_phi.set_title("phi(t)"); ax_phi.grid(True)
    ax_r   = fig1.add_subplot(gs[2,1]); ax_r.set_title("r(t)");   ax_r.grid(True)
    ax_th  = fig1.add_subplot(gs[3,0]); ax_th.set_title("theta(t)"); ax_th.grid(True)
    ax_unused= fig1.add_subplot(gs[3,1]); ax_unused.axis("off")

    time_offset = 0
    for i, phase_arr in enumerate(phases_data):
        color_ = phases_colors[i]
        n_s, dim_s = phase_arr.shape
        t_s = np.arange(n_s) + time_offset

        if dim_s>=2:
            ax_xy.plot(phase_arr[:,0], phase_arr[:,1], color=color_, label=f"Phase{i} (dim={dim_s})")
        if dim_s>0:
            ax_y.plot(t_s, phase_arr[:,0], color=color_)
        if dim_s>1:
            ax_z.plot(t_s, phase_arr[:,1], color=color_)
        if dim_s>2:
            ax_phi.plot(t_s, phase_arr[:,2], color=color_)
        if dim_s==12:
            ax_r.plot(t_s, phase_arr[:,3], color=color_)
            ax_th.plot(t_s, phase_arr[:,4], color=color_)
        elif dim_s==10:
            ax_th.plot(t_s, phase_arr[:,3], color=color_)
        time_offset += n_s

    ax_xy.legend()

    # ------------------------------------------------
    # FIGURE 2: configuration plot (modified)
    # ------------------------------------------------
    fig2, ax2 = plt.subplots(figsize=(8,6))
    if not save:
        ax2.set_title(f"{plot_title} - Drone and Pendulum Configuration")
    ax2.set_xlabel("y [m]")
    ax2.set_ylabel("z [m]")
    ax2.grid(True)
    ax2.set_ylim(-0.5, 3)  # Lower y limit increased to -0.5
    ax2.set_aspect('equal', adjustable='box')

    # --- Plot target state in orange dot ---
    if reference_xy is not None:
        # Use dot marker and orange color, with a high zorder to appear on top
        ax2.plot(reference_xy[0], reference_xy[1], marker='o',  linestyle='None', color='orange', markersize=8, label="Target State", zorder=10)
    
    # --- Plot initial state in grey if provided ---
    if initial_state is not None:
        init = np.array(initial_state).ravel()
        ax2.plot(init[0], init[1], marker='o', color='grey', markersize=10, label="Initial State", zorder=10)
    
    # --- Plot the number in the upper left corner with reduced font sizes ---
    if number is not None:
        ax2.text(0.05, 0.95, f"{number})", transform=ax2.transAxes,
                 horizontalalignment='left', verticalalignment='top',
                 fontsize=20, bbox=dict(facecolor='white', alpha=1.0, edgecolor='none', pad=4))
    
    l0 = mpc.opts.l0
    L  = mpc.opts.L_rot

    time_offset_conf = 0
    for i, phase_arr in enumerate(phases_data):
        color_ = phases_colors[i]
        n_s, dim_s = phase_arr.shape

        for idx_ in range(0, n_s, step_pose):
            row = phase_arr[idx_]
            y_   = row[0] if dim_s>0 else 0
            z_   = row[1] if dim_s>1 else 0
            phi_ = row[2] if dim_s>2 else 0

            # Draw drone body: line in black
            dy_ = L * np.cos(phi_)
            dz_ = L * np.sin(phi_)
            ax2.plot([y_ - dy_, y_ + dy_], [z_ - dz_, z_ + dz_], color='k', linewidth=2)
            # Drone center marker (red circle as before)
            ax2.plot(y_, z_, 'ro', markersize=3)
            
            # Draw pendulum line in green
            if dim_s == 12:
                r_  = row[3]
                th_ = row[4]
                load_y = y_ + r_ * np.sin(th_)
                load_z = z_ - r_ * np.cos(th_)
                ax2.plot([y_, load_y], [z_, load_z], color='green', linewidth=2)
            elif dim_s == 10:
                th_ = row[3]
                load_y = y_ + l0 * np.sin(th_)
                load_z = z_ - l0 * np.cos(th_)
                ax2.plot([y_, load_y], [z_, load_z], color='green', linewidth=2)
        time_offset_conf += n_s

    # --- Custom legend entries ---
    if legend:
        # Get current handles and labels
        handles, labels = ax2.get_legend_handles_labels()
        # Add custom handles: black line for Drone, green line for Pendulum
        drone_handle = Line2D([0], [0], color='k', lw=2, label='Drone')
        pendulum_handle = Line2D([0], [0], color='green', lw=2, label='Pendulum')
        handles.extend([drone_handle, pendulum_handle])
        legend_fontsize = 14
        ax2.legend(handles=handles, loc='lower right', fontsize=14)

    # Save configuration plot if requested
    if save:
        plot_dir = get_dir("plots")
        filename = plot_dir / f"drone/configurations_plot_{number}.pgf"
        fig2.savefig(filename, bbox_inches="tight")

    # ------------------------------------------------
    # FIGURE 3: Thrust & input signals (unchanged)
    # ------------------------------------------------
    if u_traj is not None:
        fig3 = plt.figure(figsize=(8,6))
        gs3 = fig3.add_gridspec(2,1)
        ax_thrust = fig3.add_subplot(gs3[0,0])
        ax_thrust.set_title("Left and Right Thrust vs. time")
        ax_thrust.grid(True)

        ax_inputs = fig3.add_subplot(gs3[1,0])
        ax_inputs.set_title("Control signals (dw1,dw2)")
        ax_inputs.grid(True)

        time_offset_u = 0
        for i, phase_arr in enumerate(phases_data):
            c = phases_colors[i]
            if i < len(phases_u):
                u2d = phases_u[i]
                cu  = phases_u_colors[i]
            else:
                u2d, cu = None, c

            n_s, ds = phase_arr.shape
            t_s = np.arange(n_s) + time_offset_u
            if u2d is not None:
                nu, du_ = u2d.shape
                t_u = np.arange(nu) + time_offset_u
            else:
                nu, du_ = 0, 0
                t_u = []

            c_ = getattr(mpc.opts, 'c', 1.0)

            if ds >= 12:
                w1 = phase_arr[:,10]
                w2 = phase_arr[:,11]
                thr_l = c_ * w1
                thr_r = c_ * w2
                ax_thrust.plot(t_s, thr_l, color=c, label=f"Ph{i}: left thr")
                ax_thrust.plot(t_s, thr_r, color=c, linestyle='--', label=f"Ph{i}: right thr")
            elif ds == 10:
                w1 = phase_arr[:,8]
                w2 = phase_arr[:,9]
                thr_l = c_ * w1
                thr_r = c_ * w2
                ax_thrust.plot(t_s, thr_l, color=c, label=f"Ph{i}: left thr")
                ax_thrust.plot(t_s, thr_r, color=c, linestyle='--', label=f"Ph{i}: right thr")
            else:
                if (u2d is not None) and du_ >= 2:
                    F1 = u2d[:,0]
                    F2 = u2d[:,1]
                    ax_thrust.plot(t_u, F1, color=c, label=f"Ph{i}: F1")
                    ax_thrust.plot(t_u, F2, color=c, linestyle='--', label=f"Ph{i}: F2")

            if (u2d is not None) and du_ >= 2:
                ax_inputs.plot(t_u, u2d[:,0], color=cu, label=f"Ph{i}: u0")
                ax_inputs.plot(t_u, u2d[:,1], color=cu, linestyle='--', label=f"Ph{i}: u1")

            time_offset_u += max(n_s, nu)

        ax_thrust.legend()
        ax_inputs.legend()

    # ------------------------------------------------
    # FIGURE 4: Plot of r(t) over time (new)
    # ------------------------------------------------
    if latexify:
        latexify_plot(fontsize=24)
    fig4, ax4 = plt.subplots(figsize=(8,6))
    if not save:
        ax4.set_title(f"{plot_title} - r(t) over time")
    ax4.set_xlabel("Time step")
    ax4.set_ylabel("r(t) [m]")
    ax4.grid(True)
    
    # --- Plot the number in the upper left corner for r(t) plot with reduced font size ---
    if number is not None:
        ax4.text(0.05, 0.95, f"{number})", transform=ax4.transAxes,
                 horizontalalignment='left', verticalalignment='top',
                 fontsize=32, bbox=dict(facecolor='white', alpha=1.0, edgecolor='none', pad=4))
    
    time_offset_r = 0
    for i, phase_arr in enumerate(phases_data):
        color_ = phases_colors[i]
        n_s, dim_s = phase_arr.shape
        t = np.arange(n_s) + time_offset_r
        if dim_s == 12:
            r_values = phase_arr[:, 3]
            ax4.plot(t, r_values, color=color_, label=f"Phase {i}")
        time_offset_r += n_s
    if save:
        plot_dir = get_dir("plots")
        filename = plot_dir / f"drone/r_over_time_plot_{number}.pgf"
        fig4.savefig(filename, bbox_inches="tight")

    plt.show()










