import matplotlib.pyplot as plt
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import mujoco
from trunk.trunk_utils import get_ee_position
from matplotlib.cm import tab10
import plotly.graph_objects as go
import os
import matplotlib.pyplot as plt
import numpy as np
from acados_template import latexify_plot
from utils_shared import get_dir

def plot_open_loop_plan_ee(trunkMPC, circle_x, circle_z, n_high, ee_side_id, ground_truth_mujoco=True, mjData=None, mjModel=None):
    # Extract planned x and z positions
    traj_x_ee = trunkMPC.get_planned_ee_trajectory()
    planned_x_positions = traj_x_ee[:, 0]
    planned_z_positions = traj_x_ee[:, 1]

    # Create a new plot
    plt.figure()
    plt.scatter(planned_x_positions[0], planned_z_positions[0], label='Start Point', color='red')
    plt.plot(planned_x_positions, planned_z_positions, label='Planned Trajectory', color='blue')
    plt.plot(circle_x, circle_z, '--', label='Reference Trajectory', color='green')  # Reference circle for comparison
    # Plot x-axis constraints (vertical lines)
    if trunkMPC.ub_x_ee[0] < 1e15:
        plt.axvline(x=trunkMPC.ub_x_ee[0], color='r', linestyle='--', linewidth=2, label="ub_x")
    if trunkMPC.lb_x_ee[0] > -1e15:
        plt.axvline(x=trunkMPC.lb_x_ee[0], color='r', linestyle='--', linewidth=2, label="lb_x")

    # Plot y-axis constraints (horizontal lines)
    if trunkMPC.ub_x_ee[1] < 1e15:
        plt.axhline(y=trunkMPC.ub_x_ee[1], color='b', linestyle='--', linewidth=2, label="ub_y")
    if trunkMPC.lb_x_ee[1] > -1e15:
        plt.axhline(y=trunkMPC.lb_x_ee[1], color='b', linestyle='--', linewidth=2, label="lb_y")

    if trunkMPC.has_elliptical_constraints:
        num_ellipses = trunkMPC.ellipse_centers.shape[0]
        theta = np.linspace(0, 2 * np.pi, 100)  # Parameter for drawing ellipses

        for i in range(num_ellipses):
            x_center, y_center = trunkMPC.ellipse_centers[i]
            a, b = trunkMPC.ellipse_half_axes[i]

            # Parametric equations of an ellipse
            x_ellipse = x_center + a * np.cos(theta)
            y_ellipse = y_center + b * np.sin(theta)

            plt.plot(x_ellipse, y_ellipse, 'm--', linewidth=2, label=f"Ellipsoid {i+1}" if i == 0 else None)

    if ground_truth_mujoco:

        acados_offset = -0.5

        # Deepcopy mjData and mjModel
        mjData_copy = deepcopy(mjData)
        mjModel_copy = deepcopy(mjModel)

        # Initialize ground truth trajectory
        ground_truth_x = []
        ground_truth_z = []

        steady_state_z_ee_acados_coordinates = [0.2 - acados_offset]
        mujoco.mj_forward(mjModel_copy, mjData_copy)
        # ee_pos = get_ee_position(mjData_copy, ee_side_id, steady_state_z_ee_acados_coordinates)
        # ground_truth_x.append(ee_pos[0]) 
        # ground_truth_z.append(ee_pos[2])  

        # Simulate the trajectory
        if isinstance(trunkMPC.n,list) and trunkMPC.n[-1] == 'p':
            ground_truth_steps = len(planned_x_positions[:-trunkMPC.N_list[-1]]) - trunkMPC.n_phases
        else:
            ground_truth_steps = len(planned_x_positions) - trunkMPC.n_phases
        for i in range(ground_truth_steps):
            # Record the current end-effector position
            ee_pos = get_ee_position(mjData_copy, ee_side_id, steady_state_z_ee_acados_coordinates)
            ground_truth_x.append(ee_pos[0])  
            ground_truth_z.append(ee_pos[2])  

            # Apply control and simulate forward
            control_input = trunkMPC.acados_ocp_solver.get(i, "u")

            # time step 
            time_step = trunkMPC.dt[i] 
            if control_input.size > 0:
                mjData_copy.ctrl[:] = np.array([control_input[0]] * (n_high // 2) + [control_input[1]] * (n_high // 2))
                for _ in range(int(time_step / mjModel_copy.opt.timestep)):
                    # Step simulation
                    mujoco.mj_step(mjModel_copy, mjData_copy)

        # Plot ground truth trajectory
        plt.plot(ground_truth_x, ground_truth_z, label='Ground Truth Trajectory', color='orange')
    
    # Plot details
    plt.xlabel('x position (m)')
    plt.ylabel('z position (m)')
    plt.title('Planned Trajectory vs Reference Circle')
    plt.axis('equal')  # Make the axes equal
    plt.legend()
    plt.grid()
    plt.show()

"""For generation of the pareto outlines of the sweep on the Euler cluster"""

# filter functions
def fixed_step(opts):
    return hasattr(opts, "dt") and all(x == opts.dt[0] for x in opts.dt)

def has_point_mass_model(opts):
    return hasattr(opts, "n") and "p" in opts.n

def create_pareto_frontier_outline_plot(results, filter_conditions, time_per_iter=False,
                                        show_points=False, filename="pareto.pdf", legend=False):
    """
    Creates a Pareto plot (Avg Solve Time vs. Avg Closed-Loop Cost) with optional
    scatter points and frontier outlines for each filter group. Frontier lines
    are extended up and right to fixed axis limits.

    Args:
        results (list): List of result dicts with keys: 'costs', 'solve_times',
                        'SQP_iters', 'options', 'constraint_violations'.
        filter_conditions (list): List of tuples (condition_func, color, label)
                                  defining the filtering and coloring.
        time_per_iter (bool): Use average solve time per SQP iteration if True.
        show_points (bool): If True, scatter points are plotted. If False, only
                            the frontiers and extensions are shown.
    """
    X_LIM = 1.0
    Y_LIM = 0.01

    # Collect points for each filter condition
    filtered_points = {label: [] for (_, _, label) in filter_conditions}

    for result in results:
        costs = result.get('costs', [])
        solve_times = result.get('solve_times', [])
        SQP_iters = result.get('SQP_iters', [])
        options = result.get('options', None)
        constraint_violations = result.get('constraint_violations', 0)

        avg_cost = np.mean(costs) if costs else float('inf')

        # Skip invalid points for which the solver failed
        if (options is None or len(solve_times) == 0 or
            constraint_violations > 0.01):
            continue

        # Solve time per iteration (if applicable)
        if time_per_iter:
            iters_array = np.array(SQP_iters, dtype=float)
            iters_array[iters_array == 0] = 1e-8
            solve_time_array = np.array(solve_times, dtype=float) / iters_array
        else:
            solve_time_array = np.array(solve_times, dtype=float)

        avg_solve_time = np.mean(solve_time_array)

        # Add to matching group (can belong to multiple)
        for condition_func, color, label in filter_conditions:
            if condition_func(options):
                filtered_points[label].append((avg_solve_time, avg_cost, color))

    # Helper to compute Pareto frontier (non-dominated, lower-left)
    def compute_pareto_frontier(xy_list):
        xy_sorted = sorted(xy_list, key=lambda p: p[0])
        frontier = []
        min_y = float('inf')
        for x, y, _ in xy_sorted:
            if y < min_y:
                frontier.append((x, y))
                min_y = y
        return frontier

    # Begin plot
    latexify_plot(fontsize=18)
    plt.figure(figsize=(12.2, 4.5))

    for (_, color, label) in filter_conditions:
        pts = filtered_points[label]
        if not pts:
            continue

        if show_points:
            x_vals = [p[0] for p in pts]
            y_vals = [p[1] for p in pts]
            plt.scatter(x_vals, y_vals, color=color, alpha=0.5, marker='o')

        # Compute and plot Pareto frontier
        frontier = compute_pareto_frontier(pts)
        if frontier:
            fx = [p[0] for p in frontier]
            fy = [p[1] for p in frontier]
            plt.plot(fx, fy, color=color, linewidth=2, label=label)

            # Vertical extension from leftmost point
            leftmost_x, leftmost_y = frontier[0]
            plt.plot([leftmost_x, leftmost_x], [leftmost_y, Y_LIM], color=color, linewidth=2)

            # Horizontal extension from bottommost point
            bottommost_x, bottommost_y = min(frontier, key=lambda p: p[1])
            plt.plot([bottommost_x, X_LIM], [bottommost_y, bottommost_y], color=color, linewidth=2)

    # Final styling
    #plt.title('Pareto Plot with Frontier Outlines')
    plt.xlabel('Average Solve Time (s)')
    plt.ylabel('Average Closed-Loop Cost')
    plt.xlim(0, 0.09)
    plt.ylim(0.007, 0.0076)
    plt.grid(True)
    if legend:
        plt.legend(
                loc='upper right',
                labelspacing=0.2,
                handlelength=1.5,
                handletextpad=0.4,
                borderaxespad=0.3,
                frameon=True,
                fontsize=18
            )
    #plt.tight_layout()
    trunk_plot_dir = get_dir("plots/trunk")
    filename = trunk_plot_dir / filename
    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight')
    plt.show()


""" Less relevant plotting functions"""

def create_pareto_plots(results, time_per_iter=True):
    """
    Creates two separate plots:
    - Pareto plot of average solve times vs. average closed-loop costs (with error bars).
    - Pareto plot of constraint violations vs. solve times.
    Args:
        results: List of dictionaries containing 'costs', 'solve_times', 'constraint_violations', and 'options'.
    """
    avg_costs = []
    constraint_violations = []
    avg_solve_times = []
    std_solve_times = []
    labels = []

    # Extract average costs, solve times, standard deviations, and labels
    for result in results:
        avg_cost = np.mean(result['costs'])
        solve_times = result['solve_times']
        SQP_iters = result['SQP_iters']
        if time_per_iter:
            solve_time = np.array(solve_times) / np.array(SQP_iters)
        else:
            solve_time = np.array(solve_times)
        avg_solve_time = np.mean(solve_time)
        std_solve_time = np.std(solve_time)
        constraint_violation = np.mean(result['constraint_violations'])
        options = result['options']
        
        label = f"n={options.n}, N_list={options.N_list}, solver={options.nlp_solver_type}"
        avg_costs.append(avg_cost)
        avg_solve_times.append(avg_solve_time)
        std_solve_times.append(std_solve_time)
        constraint_violations.append(constraint_violation)
        labels.append(label)

    # First plot: Average solve time vs. average closed-loop cost
    plt.figure(figsize=(8, 6))
    for i, (cost, solve_time, std_dev, label, options) in enumerate(zip(
            avg_costs, avg_solve_times, std_solve_times, labels, [r['options'] for r in results])):
        marker = 'o' if options.nlp_solver_type == "SQP" else 'D'  # 'o' for SQP, 'D' (diamond) for SQP_RTI
        
        # Scatter plot with error bars
        plt.errorbar(solve_time, cost, xerr=std_dev, fmt=marker, label=label, alpha=0.7, capsize=5)

    plt.title('Pareto Plot: Avg Solve Times vs. Avg Closed-Loop Costs')
    plt.xlabel('Average Solve Time (s)')
    plt.ylabel('Average Closed-Loop Cost')
    plt.legend(loc='best', fontsize='small', title='Simulation Parameters')
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    # Second plot: Constraint violations vs. solve time
    plt.figure(figsize=(8, 6))
    for i, (violation, solve_time, label, options) in enumerate(zip(
            constraint_violations, avg_solve_times, labels, [r['options'] for r in results])):
        marker = 'o' if options.nlp_solver_type == "SQP" else 'D'  # 'o' for SQP, 'D' (diamond) for SQP_RTI
        
        # Scatter plot
        plt.scatter(solve_time, violation, marker=marker, label=label, alpha=0.7)
    
    plt.title('Pareto Plot: Constraint Violations vs. Solve Times')
    plt.xlabel('Average Solve Time (s)')
    plt.ylabel('Constraint Violations')
    plt.legend(loc='best', fontsize='small', title='Simulation Parameters')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    from matplotlib.cm import tab10

def compute_number_of_decision_variables(options):
    total_vars = 0
    # n[i] can be int or 'p'
    # N_list[i] is the horizon length for that sub-phase
    # direct multiple shooting formula: (N_i + 1)*state_dim + N_i*input_dim

    for i in range(len(options.N_list)):
        N_i = options.N_list[i]
        # dimension
        if i < len(options.n):
            n_i = options.n[i]
            if isinstance(n_i, int):
                state_dim = 2 * n_i
            elif n_i == 'p':
                # point mass with 4 states
                state_dim = 4
            else:
                # unknown
                state_dim = 0
        else:
            state_dim = 0
        input_dim = 2
        total_vars += (N_i + 1)*state_dim + (N_i)*input_dim
    return total_vars


def create_pareto_plots_clustered(results, plot_sqp=True, plot_error_bars=True, time_per_iter=True,
                                  plot_group1=True, plot_group2=True, plot_group3=True, plot_group4=True):
    """
    Creates Pareto plots with results clustered based on model scheduling and step size rates.
    - Different groups are assigned different colors.
    - Circles (o) for SQP and Diamonds (D) for SQP_RTI solvers.
    - Plots:
        1) Average solve time vs. average closed-loop cost
        2) Constraint violation vs. average solve time
        3) Average closed-loop cost vs. number of decision variables

    Legend clarifies:
        - Group 1: single model, no step size increase
        - Group 2: single model, step size increase
        - Group 3: multi-model, no step size increase
        - Group 4: multi-model, step size increase

    Args:
        results: List of dictionaries containing 'costs', 'solve_times', 'constraint_violations', and 'options'.
        plot_sqp (bool): If True, only plot results for SQP solver; if False, only plot SQP_RTI.
        plot_error_bars (bool): If True, include error bars (solve time) on the first plot; if False, only scatter.
        time_per_iter (bool): If True, average solve time is divided by # of SQP iterations.
        plot_group1 (bool): Toggle for plotting group 1 (baseline, single model, no step size increase).
        plot_group2 (bool): Toggle for plotting group 2 (single model, increasing step size).
        plot_group3 (bool): Toggle for plotting group 3 (multi-model, no step size increase).
        plot_group4 (bool): Toggle for plotting group 4 (multi-model, increasing step size).
    """

    # We'll maintain colors and textual descriptions separately
    group_colors = {
        'group1': 'blue',
        'group2': 'green',
        'group3': 'red',
        'group4': 'purple'
    }
    group_labels = {
        'group1': 'Group 1: single model, no step size increase',
        'group2': 'Group 2: single model, step size increase',
        'group3': 'Group 3: multi-model, no step size increase',
        'group4': 'Group 4: multi-model, step size increase'
    }

    avg_costs = []
    avg_solve_times = []
    std_solve_times = []
    constraint_violations = []
    solver_types = []
    groups = []
    num_vars_list = []

    def classify_group(n_schedule, dt_array):
        rate = 1.0
        if len(dt_array) > 1 and dt_array[0] != 0:
            rate = dt_array[1] / dt_array[0]
        if len(n_schedule) == 1 and rate == 1.0:
            return 'group1'
        elif len(n_schedule) == 1 and rate > 1.0:
            return 'group2'
        elif len(n_schedule) > 1 and rate == 1.0:
            return 'group3'
        else:
            return 'group4'

    # === Extract metrics ===
    for result in results:
        print(result['costs'])
        avg_cost = np.mean(result['costs'])
        solve_times = result['solve_times']
        SQP_iters = result['SQP_iters']
        constraint_violation = np.mean(result['constraint_violations'])
        options = result['options']

        # If time_per_iter is True, divide solve time by the number of SQP iterations
        if time_per_iter:
            iters_array = np.array(SQP_iters, dtype=float)
            iters_array[iters_array == 0] = 1e-8  # avoid divide-by-zero
            solve_time_array = np.array(solve_times, dtype=float) / iters_array
        else:
            solve_time_array = np.array(solve_times, dtype=float)

        avg_solve_time = np.mean(solve_time_array)
        std_solve_time = np.std(solve_time_array)

        # Classify group
        group = classify_group(options.n, options.dt)

        # Filter out results that don't match the chosen solver
        if plot_sqp and options.nlp_solver_type != "SQP":
            continue
        if not plot_sqp and options.nlp_solver_type != "SQP_RTI":
            continue

        # Filter groups according to user-provided flags
        if group == 'group1' and not plot_group1:
            continue
        if group == 'group2' and not plot_group2:
            continue
        if group == 'group3' and not plot_group3:
            continue
        if group == 'group4' and not plot_group4:
            continue

        # Compute number of decision variables
        total_vars = compute_number_of_decision_variables(options)

        # Store data
        avg_costs.append(avg_cost)
        avg_solve_times.append(avg_solve_time)
        std_solve_times.append(std_solve_time)
        constraint_violations.append(constraint_violation)
        solver_types.append(options.nlp_solver_type)
        groups.append(group)
        num_vars_list.append(total_vars)

    solver_label = "SQP" if plot_sqp else "SQP_RTI"

    # === 1) Plot: Avg Solve Time vs. Avg Closed-Loop Cost ===
    plt.figure(figsize=(8, 6))
    for cost, solve_time, std_dev, solver, group in zip(
            avg_costs, avg_solve_times, std_solve_times, solver_types, groups):
        marker = 'o' if solver == "SQP" else 'D'
        color = group_colors[group]
        if plot_error_bars:
            plt.errorbar(solve_time, cost, xerr=std_dev, fmt=marker, alpha=0.7, capsize=5, color=color)
        else:
            plt.scatter(solve_time, cost, marker=marker, alpha=0.7, color=color)

    # Legend for groups with textual meaning
    for grp_key in group_colors:
        plt.scatter([], [], color=group_colors[grp_key], label=group_labels[grp_key], alpha=0.7)

    plt.title(f'Pareto Plot ({solver_label}): Avg Solve Times vs. Avg Closed-Loop Costs')
    plt.xlabel('Average Solve Time (s)')
    plt.ylabel('Average Closed-Loop Cost')
    plt.legend(loc='best', fontsize='small', title='Groups')
    plt.grid(True)
    plt.tight_layout()
    plt.xlim(0, 0.1)  # Set x-axis limits
    plt.ylim(0, 10) # Set y-axis limits
    plt.show()

    # === 2) Plot: Constraint Violations vs. Solve Time ===
    plt.figure(figsize=(8, 6))
    for violation, solve_time, solver, group in zip(
            constraint_violations, avg_solve_times, solver_types, groups):
        marker = 'o' if solver == "SQP" else 'D'
        color = group_colors[group]
        plt.scatter(solve_time, violation, marker=marker, alpha=0.7, color=color)

    # Legend for groups with textual meaning
    for grp_key in group_colors:
        plt.scatter([], [], color=group_colors[grp_key], label=group_labels[grp_key], alpha=0.7)

    plt.title(f'Pareto Plot ({solver_label}): Constraint Violations vs. Solve Times')
    plt.xlabel('Average Solve Time (s)')
    plt.ylabel('Constraint Violations')
    plt.legend(loc='best', fontsize='small', title='Groups')
    plt.grid(True)
    plt.tight_layout()
    plt.xlim(0, 0.1)  # Set x-axis limits
    plt.ylim(0, 10) # Set y-axis limits
    plt.show()

    # === 3) Plot: Avg Closed-Loop Cost vs. Number of Decision Variables ===
    plt.figure(figsize=(8, 6))
    for cost, num_vars, solver, group in zip(avg_costs, num_vars_list, solver_types, groups):
        marker = 'o' if solver == "SQP" else 'D'
        color = group_colors[group]
        plt.scatter(num_vars, cost, marker=marker, alpha=0.7, color=color)

    # Legend for groups with textual meaning
    for grp_key in group_colors:
        plt.scatter([], [], color=group_colors[grp_key], label=group_labels[grp_key], alpha=0.7)

    plt.title(f'Pareto Plot ({solver_label}): Avg Cost vs. Number of Decision Variables')
    plt.xlabel('Number of Decision Variables (approx)')
    plt.ylabel('Average Closed-Loop Cost')
    plt.legend(loc='best', fontsize='small', title='Groups')
    plt.grid(True)
    plt.tight_layout()
    plt.xlim(0, 1e4)  # Set x-axis limits
    plt.ylim(0, 10) # Set y-axis limits
    plt.show()

def create_interactive_cost_plots(results):
    """
    Create interactive plots showing costs over time and accumulated costs over time.
    Args:
        results (list): List of result dictionaries containing 'costs' and 'options'.
    """
    # Data lists
    costs_over_time = []
    accumulated_costs = []
    labels = []
    n_values = []

    # Process results
    for result in results:
        costs = result['costs']
        accumulated = np.cumsum(costs)
        options = result['options']
        label = f"n={options.n}, N_list={options.N_list}, solver={options.nlp_solver_type}"
        
        costs_over_time.append(costs)
        accumulated_costs.append(accumulated)
        labels.append(label)
        n_values.append(tuple(map(str, options.n)))  # Convert n to strings for uniformity

    # Unique color mapping for 'n'
    unique_n_values = sorted(set(n_values))
    color_map = {n: f"rgba({tab10(i % 10)[0]*255}, {tab10(i % 10)[1]*255}, {tab10(i % 10)[2]*255}, 0.7)" 
                 for i, n in enumerate(unique_n_values)}

    # Plot costs over time
    fig1 = go.Figure()
    for i, (costs, label, n) in enumerate(zip(costs_over_time, labels, n_values)):
        color = color_map[n]
        
        fig1.add_trace(go.Scatter(
            x=list(range(len(costs))),
            y=costs,
            mode='lines',
            line=dict(color=color),
            name=label,
            hoverinfo='name+y'
        ))

    fig1.update_layout(
        title='Costs Over Time',
        xaxis_title='Time Steps',
        yaxis_title='Cost',
        legend_title='Simulation Parameters',
        template="plotly",
        hovermode='closest'
    )

    fig1.show()

    # Plot accumulated costs over time
    fig2 = go.Figure()
    for i, (acc_costs, label, n) in enumerate(zip(accumulated_costs, labels, n_values)):
        color = color_map[n]
        
        fig2.add_trace(go.Scatter(
            x=list(range(len(acc_costs))),
            y=acc_costs,
            mode='lines',
            line=dict(color=color),
            name=label,
            hoverinfo='name+y'
        ))

    fig2.update_layout(
        title='Accumulated Costs Over Time',
        xaxis_title='Time Steps',
        yaxis_title='Accumulated Cost',
        legend_title='Simulation Parameters',
        template="plotly",
        hovermode='closest'
    )

    fig2.show()
