"""
Closed-loop simulation functionality
"""
from trunk.trunk_mpc_2d_acados import TrunkMPC2DOptions, TrunkMPC2DEETracking
from trunk.Trajectory_finding_script import select_periodic_trajectory
import mujoco
import numpy as np
from trunk.trunk_utils import compute_q
from trunk.trunk_utils import get_ee_position
from trunk.plotting_utils_trajectory_tracking import plot_open_loop_plan_ee
from trunk.Trajectory_finding_script import select_periodic_trajectory
import os
import mediapy as media
import casadi as ca
from trunk.trunk_utils import compute_q_casadi
from trunk.ODE_chains.ODE_utils import get_ode_params_dict
import pickle
import gc
from trunk.plotting_utils_trajectory_tracking import create_pareto_plots_clustered, create_interactive_cost_plots
from copy import deepcopy

def closed_loop_ee_tracking_mujoco(n, n_high, Mjmodel, Mjdata, duration, framerate, create_video, plot_open_loop_plan_bool, trajectory_x, trajectory_z, trunkMPC, controller_interval=None, perturbation_variance=0.0):
    # functions for mapping the states between chains of different dimensions, in case we try to control a higher dimensional chain with a lower dimensional model, otherwhise the mappings are the identity
    q_high_dim = ca.SX.sym(f'q_high_dim_{n_high}', n_high)  # Symbolic q_high for 16-link
    q_dot_high_dim = ca.SX.sym(f'q_dot_high_dim_{n_high}', n_high)  # Symbolic q_high_dot for 16-link

    ode_params_dict = get_ode_params_dict()
    link_length_high_dim = 2*ode_params_dict[f"{n_high}"]["l"]
    link_length_low_dim = 2*ode_params_dict[f"{n[0]}"]["l"]
    q_low_dim = compute_q_casadi(q_high_dim, link_length_high_dim, link_length_low_dim, n[0])
    jacobian_q_low_q_high = ca.jacobian(q_low_dim, q_high_dim)  # Jacobian of q_8 w.r.t. q_16
    q_dot_low_dim = jacobian_q_low_q_high @ q_dot_high_dim  # Velocity mapping for 8-link
    q_dot_low_dim_fn = ca.Function('q_dot_low_dim_fn', [q_high_dim, q_dot_high_dim], [q_dot_low_dim])  # CasADi function

    # set control frequency to the timestep of the MPC
    if controller_interval is not None:
        controller_interval = controller_interval
    else:
        controller_interval = trunkMPC.dt
    
    # simulate
    ee_side_id = mujoco.mj_name2id(Mjmodel, mujoco.mjtObj.mjOBJ_SITE, "endeffector")
    frames = []
    control_inputs = []
    ee_pos = []
    costs = []
    solve_times = []
    SQP_iters = []
    constraint_violations = 0
    last_update_time = -controller_interval  # Initialize to ensure the controller is applied at the start
    with mujoco.Renderer(Mjmodel, width=1920, height=1080) as renderer:
        while Mjdata.time < duration:
            # extract and save endeffector position
            mujoco.mj_forward(Mjmodel, Mjdata) # just to be save that all positions etc are updated
            #print(f"MuJoCo Timestep: {Mjmodel.opt.timestep}")

            ee_pos.append(get_ee_position(Mjdata, ee_side_id, [0.7]))

            x_ee, _, y_ee = ee_pos[-1]  # Assuming ee_pos is a list of (x, y) positions

            # Check violations for x position
            if x_ee > trunkMPC.ub_x_ee[0]:  
                constraint_violations += (x_ee - trunkMPC.ub_x_ee[0])
            if x_ee < trunkMPC.lb_x_ee[0]:  
                constraint_violations += (trunkMPC.lb_x_ee[0] - x_ee)

            # Check violations for y position
            if y_ee > trunkMPC.ub_x_ee[1]:  
                constraint_violations += (y_ee - trunkMPC.ub_x_ee[1])
            if y_ee < trunkMPC.lb_x_ee[1]:  
                constraint_violations += (trunkMPC.lb_x_ee[1] - y_ee)

            # Apply controller at given frequency
            if Mjdata.time - last_update_time >= controller_interval:
                x0 = np.concatenate((Mjdata.qpos, Mjdata.qvel))
                q0_high_dim = Mjdata.qpos
                q0_dot_high_dim = Mjdata.qvel
                q0_mapped = compute_q(q0_high_dim, link_length_high_dim, link_length_low_dim, n[0])
                q0_dot_mapped = q_dot_low_dim_fn(q0_high_dim, q0_dot_high_dim).full().flatten()
                x0_mapped = np.concatenate((q0_mapped, q0_dot_mapped))
                trunkMPC.set_initial_state(x0_mapped)
                trunkMPC.set_reference_trajectory(Mjdata.time)
                u0, solve_time_tot, sqp_iters = trunkMPC.solve_ocp()
                step = int(round(Mjdata.time / trunkMPC.dt[0]))
                if step % 200 == 0 and plot_open_loop_plan_bool:
                    plot_open_loop_plan_ee(trunkMPC, trajectory_x, trajectory_z, n_high, ee_side_id, ground_truth_mujoco=True, mjData=Mjdata, mjModel=Mjmodel)
                u0_mapped = u0#compute_controls_inverse(u0[0], u0[1], link_length_high_dim, link_length_low_dim, ode_params_dict[f"{n_high}"]["c"], ode_params_dict[f"{n}"]["c"], 2, n_high, n)#u0
                Mjdata.ctrl[:] = np.array([u0_mapped[0]] * (n_high // 2) + [u0_mapped[1]] * (n_high // 2))
                last_update_time = Mjdata.time

                # append solve time
                solve_times.append(solve_time_tot)
                SQP_iters.append(sqp_iters)
            
            # compute costs
            if trunkMPC.multi_phase:
                Q, R = trunkMPC.ocp.cost[0].Q, trunkMPC.ocp.cost[0].R
            else:
                Q, R = trunkMPC.ocp.cost.Q, trunkMPC.ocp.cost.R
            ee_pos_diff = np.delete(ee_pos[-1], 1) - trunkMPC.Phi_t(Mjdata.time)
            site_velocity = np.zeros(6)
            mujoco.mj_objectVelocity(Mjmodel, Mjdata, mujoco.mjtObj.mjOBJ_SITE, ee_side_id, site_velocity, flg_local=False)
            ee_site_linear_velocity = site_velocity[3:]
            ee_dot_diff = np.delete(ee_site_linear_velocity,1) - trunkMPC.Phi_dot_t(Mjdata.time)
            y_diff = np.concatenate((ee_pos_diff, ee_dot_diff))
            costs.append(y_diff.T @ Q @ y_diff + u0_mapped.T @ R @  u0_mapped)

            mujoco.mj_step(Mjmodel, Mjdata)
            # add noise to the state 
            # Add random noise to joint velocities (or positions, if needed)
            noise = np.random.normal(0.000, np.sqrt(perturbation_variance), Mjdata.qpos.shape)
            # Mjdata.qvel += noise  # Apply noise
            # Mjdata.qpos += noise  # Apply noise to positions
            # try if applying a force at the endeffector produces more predictable results
            endeffector_body_id = Mjmodel.site_bodyid[ee_side_id]
            sigma = np.sqrt(perturbation_variance)
            force_2d = np.random.normal(0, sigma, size=2)
            force_x, force_z = force_2d
            Mjdata.xfrc_applied[endeffector_body_id, 0] = force_x
            Mjdata.xfrc_applied[endeffector_body_id, 1] = 0.0
            Mjdata.xfrc_applied[endeffector_body_id, 2] = force_z

            control_inputs.append(u0_mapped)

            if create_video and len(frames) < Mjdata.time * framerate:
                renderer.update_scene(Mjdata, camera="yz_view")
                pixels = renderer.render()
                frames.append(pixels)

    if create_video:
        video_folder = os.path.join(os.path.dirname(__file__), "trajectories", "videos", f"reference_tracking_2d_trunk_{n_high}_{n}.mp4")
        media.write_video(video_folder, frames, fps=framerate)
        print('Rendering done')

    return control_inputs, ee_pos, costs, solve_times, SQP_iters, constraint_violations


def closed_loop_ee_tracking_acados(n, n_high, duration, plot_open_loop_plan_bool, trajectory_x, trajectory_z, trunkMPC, trunkMPC_switching=False, switching_time=None):  
    "closed loop sim in acados"
    
    control_inputs = []
    ee_pos = []
    costs = []
    solve_times = []
    SQP_iters = []
    constraint_violations = 0
    time = 0
    
    x0 = np.zeros((2 * n[0],))

    states = []
    time_list = []

    # Prepare EE-velocity function 
    q = ca.MX.sym('q', n[0])  # generalized coodinates
    q_dot = ca.MX.sym('q_dot', n[0])  # generalized coodinates
    
    # map high dimensional q and q_dots to endeffector position
    x_pos_ee, y_pos_ee = trunkMPC.forward_kinematics_casadi(q, n[0])   
    x_ee_dot = ca.jacobian(x_pos_ee, q) @ q_dot 
    y_ee_dot = ca.jacobian(y_pos_ee, q) @ q_dot   

    EE_dot_fct = ca.Function('ee_dot', [q, q_dot], [ca.vertcat(x_ee_dot, y_ee_dot)])

    if switching_time is not None:
        # save mapping between states
        n_high_dim = trunkMPC_switching.n[0]
        n_low_dim = trunkMPC_switching.n[1]
        q_high_dim = ca.MX.sym('q', n_high_dim)  # generalized coordinates
        q_dot_high_dim = ca.MX.sym('q_dot', n_high_dim)  # generalized velocities
        
        # map high dimensional q and q_dots to low dimensional ones
        link_length_high_dim = 2 * trunkMPC_switching.ode_params_dict[f"{n_high_dim}"]["l"]
        link_length_low_dim = 2 * trunkMPC_switching.ode_params_dict[f"{n_low_dim}"]["l"]
        q_low_dim = compute_q_casadi(q_high_dim, link_length_high_dim, link_length_low_dim, n_low_dim)
        jacobian_q_low_dim_q_high_dim = ca.jacobian(q_low_dim, q_high_dim)
        q_dot_low_dim = jacobian_q_low_dim_q_high_dim @ q_dot_high_dim  # Velocity mapping for high-dimensional system
        x_low_dim_fct = ca.Function('x_low_dim_fct', [ca.vertcat(q_high_dim, q_dot_high_dim)], [ca.vertcat(q_low_dim, q_dot_low_dim)])  # CasADi function

        acados_sim_switched = trunkMPC_switching.create_sim_solver(trunkMPC_switching.multi_phase, stage=2)


    acados_sim = trunkMPC.acados_sim_solver
    switched = False
    while time < duration:
    
        trunkMPC.set_reference_trajectory(time)
        trunkMPC.set_initial_state(x0)
        u0, solve_time_tot, sqp_iters = trunkMPC.solve_ocp()

        ee_pos.append(trunkMPC.forward_kinematics_casadi(x0[:n[0]], n[0]))
        # Extract end-effector position
        x_ee, y_ee = ee_pos[-1]  # Assuming ee_pos is a list of (x, y) positions

        # Check violations for x position
        if x_ee > trunkMPC.ub_x_ee[0]:  
            constraint_violations += (x_ee - trunkMPC.ub_x_ee[0])
        if x_ee < trunkMPC.lb_x_ee[0]:  
            constraint_violations += (trunkMPC.lb_x_ee[0] - x_ee)

        # Check violations for y position
        if y_ee > trunkMPC.ub_x_ee[1]:  
            constraint_violations += (y_ee - trunkMPC.ub_x_ee[1])
        if y_ee < trunkMPC.lb_x_ee[1]:  
            constraint_violations += (trunkMPC.lb_x_ee[1] - y_ee)
        
        acados_sim.set("x", x0)
        acados_sim.set("u", u0)
        status = acados_sim.solve()
        if status != 0:
            raise RuntimeError(f"ACADOS sim solver failed with status {status}.")
        
        step = int(round(time / trunkMPC.dt[0]))
        if step % 250 == 0 and plot_open_loop_plan_bool:
            print("computation time:", solve_time_tot)
            plot_open_loop_plan_ee(trunkMPC, trajectory_x, trajectory_z, n_high, None, ground_truth_mujoco=False)

        # compute costs
        if trunkMPC.multi_phase:
            Q, R = trunkMPC.ocp.cost[0].Q, trunkMPC.ocp.cost[0].R
        else:
            Q, R = trunkMPC.ocp.cost.Q, trunkMPC.ocp.cost.R
        ee_pos_diff = np.array(ee_pos[-1]) - trunkMPC.Phi_t(time)
        ee_linear_velocity = EE_dot_fct(x0[:n[0]], x0[-n[0]:]).full().squeeze()
        #test = trunkMPC.Phi_dot_t(time)
        ee_dot_diff = ee_linear_velocity - trunkMPC.Phi_dot_t(time)
        y_diff = np.concatenate((ee_pos_diff, ee_dot_diff))
        costs.append(y_diff.T @ Q @ y_diff + u0.T @ R @  u0)
        
        x_next = acados_sim.get("x")
        states.append(x0)
        x0 = x_next
        time_list.append(time)
        time += trunkMPC.dt[0]

        control_inputs.append(u0)
        solve_times.append(solve_time_tot)
        SQP_iters.append(sqp_iters)

        # switch trunk MPC after given time
        if switching_time is not None and not switched and time > switching_time:
            # try better initialization, this is currently not very generla, just for testing
            # get previous trajetory of MPC with shorter horizon but high fidelity model
            traj_x_prev, traj_u_prev = trunkMPC.get_planned_trajectory()
            traj_x_guess = traj_x_prev[:trunkMPC_switching.N_list[0]+1]
            traj_u_guess = traj_u_prev[:trunkMPC_switching.N_list[0]]
            x_1_end = traj_x_prev[trunkMPC_switching.N_list[0]] # x-rest should include x_last and all remaining states by simulating zeros forward
            link_length_high_dim = 2*trunkMPC_switching.ode_params_dict[f"{trunkMPC.n}"]["l"]
            link_length_low_dim = 2*trunkMPC_switching.ode_params_dict[f"{trunkMPC_switching.n[1]}"]["l"]
            x_low_dim = x_low_dim_fct(x_1_end)
            traj_x_guess.append(x_low_dim.full().squeeze())
            traj_u_guess.append(np.empty((0, 0)))
            for i in range(trunkMPC_switching.N_list[1]):
                acados_sim_switched.set("x", traj_x_guess[-1])
                acados_sim_switched.set("u", np.zeros(2))
                status = acados_sim_switched.solve()
                if status != 0:
                    raise RuntimeError(f"ACADOS low dimensional sim solver failed with status {status}.")
                traj_x_guess.append(acados_sim_switched.get("x"))
                traj_u_guess.append(np.zeros(2))

            # map x_rest to lower dimensional state after trunkMPC.N_list[0] 
            if True:
                trunkMPC_switching.set_initial_guess(traj_x_guess, traj_u_guess)
            trunkMPC = trunkMPC_switching
            switched = True
            print("Switch Trunk MPC")

    return control_inputs, states, time_list, ee_pos, costs, solve_times, SQP_iters, constraint_violations
    
