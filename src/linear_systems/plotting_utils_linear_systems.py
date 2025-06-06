#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, QhullError
from acados_template import latexify_plot
import control as ct
from utils_shared import get_dir

def create_filtered_trace_increase_plot(results, filter_conditions,
                                       default_color='lightgray', default_label='Other Results',
                                       plot_uncondensed=True, plot_condensed=False,
                                       pareto_frontier=True,
                                       frontier_linewidth=2.0,
                                       point_alpha=0.3, filename="trace_pareto"):
    """
    Plots percentage increase in trace(P) vs. number of optimization variables.

    * If pareto_frontier=True:
      - We compute a "bottom frontier" of the convex hull for each group (uncondensed/condensed).
      - We plot those hull lines in a stepped manner (the "lowest boundary").
      - We also add two dashed lines (horizontal & vertical) from the minimum-y and minimum-x points.
      - Additionally, we scatter the points in the background (with transparency 'point_alpha').
    * If pareto_frontier=False, we produce the original scatter plot only.
    
    Additional adjustable parameters:
      - frontier_linewidth: thickness of the Pareto frontier lines
      - point_alpha: transparency of points when the Pareto frontier is drawn

    Points with y < 0 are discarded.
    """
    filtered_points = {label: {"uncondensed": [], "condensed": []}
                       for (_, _, label) in filter_conditions}
    default_points = {"uncondensed": [], "condensed": []}

    # ------------------------------------------------
    # 1) Gather & filter the (x, y) data by label
    # ------------------------------------------------
    for result in results:
        P_cont = result.get('P_cont')
        opts = result.get('opts')

        m = opts["m"]
        dt = opts["dt"]
        A = result.get('A')
        B = result.get('B')

        # Hardcoded Q, R
        Q = np.eye(A.shape[0])
        R = 0.1 * np.eye(m)

        # Discretize system
        sys_gt_ct = ct.ss(A, B, np.eye(A.shape[0]), np.zeros((A.shape[0], B.shape[1])))
        sys_gt_disc_step = ct.c2d(sys_gt_ct, dt, method='zoh')
        A_d = sys_gt_disc_step.A
        B_d = sys_gt_disc_step.B

        P_LQR, _, _ = ct.dare(A_d, B_d, Q * dt, R * dt)
        if P_LQR is None:
            continue

        trace_P_cont = np.trace(P_cont)
        trace_P_LQR = np.trace(P_LQR)
        perc_increase = (trace_P_cont - trace_P_LQR) / trace_P_LQR * 100.0
        # Skip negative y-values, arise due to numerical issues, like due to instability
        if perc_increase < 0:
            continue

        x_uncond = result.get('N_optim_vars_uncondensed', None)
        x_cond = result.get('N_optim_vars_condensed', None)

        # Skip if we can't plot it
        if (plot_uncondensed and x_uncond is None) or (plot_condensed and x_cond is None):
            continue

        # Decide which label-group it belongs to
        matched_any = False
        for condition_func, color, label in filter_conditions:
            if condition_func(result):
                if plot_uncondensed:
                    filtered_points[label]["uncondensed"].append((x_uncond, perc_increase))
                if plot_condensed:
                    filtered_points[label]["condensed"].append((x_cond, perc_increase))
                matched_any = True
                break
        if not matched_any:
            if plot_uncondensed:
                default_points["uncondensed"].append((x_uncond, perc_increase))
            if plot_condensed:
                default_points["condensed"].append((x_cond, perc_increase))

    # ------------------------------------------------
    # 2) Function to find & plot the "bottom frontier"
    #    plus horizontal & vertical lines
    # ------------------------------------------------
    def plot_bottom_frontier(pts, color, lbl):
        """
        1) Compute the hull for pts,
        2) Extract the lower boundary from leftmost to rightmost in hull CCW order,
        3) Plot it as a single line,
        4) Then from the boundary's min-y point, draw a horizontal line to max x,
        5) Then from the boundary's min-x point, draw a vertical line to max y.
        """
        if len(pts) < 3:
            # If fewer than 3 points, we can't form a hull.
            return

        arr = np.array(pts)  # shape: (N,2)
        try:
            hull = ConvexHull(arr)
        except QhullError:
            return

        hull_indices = hull.vertices  # these are indices into arr (CCW order)
        xs_hull = arr[hull_indices, 0]

        # Identify leftmost hull vertex (lowest x) and rightmost hull vertex (highest x)
        leftmost_in_hullarray = np.argmin(xs_hull)
        rightmost_in_hullarray = np.argmax(xs_hull)

        # Walk the hull from leftmost to rightmost
        def walk_hull(start_idx, end_idx):
            chain = []
            i = start_idx
            while True:
                real_index = hull_indices[i]
                chain.append(real_index)
                if i == end_idx:
                    break
                i = (i + 1) % len(hull_indices)
            return chain

        forward_chain = walk_hull(leftmost_in_hullarray, rightmost_in_hullarray)
        backward_chain = walk_hull(rightmost_in_hullarray, leftmost_in_hullarray)

        forward_pts = arr[forward_chain, :]
        backward_pts = arr[backward_chain, :]

        # Pick the chain with the smaller average y => the "bottom" chain
        if np.mean(forward_pts[:, 1]) < np.mean(backward_pts[:, 1]):
            bottom_chain_pts = forward_pts
        else:
            bottom_chain_pts = backward_pts

        # Plot the lower boundary as a single line
        plt.plot(bottom_chain_pts[:, 0], bottom_chain_pts[:, 1],
                 color=color, linewidth=frontier_linewidth, label=lbl)

        # Identify min-y, min-x from that boundary, plus max-x, max-y
        x_vals = bottom_chain_pts[:, 0]
        y_vals = bottom_chain_pts[:, 1]

        x_min = np.min(x_vals)
        x_max = np.max(x_vals)
        y_min = np.min(y_vals)
        y_max = np.max(y_vals)

        # Find point with min-y
        bottom_mask = np.isclose(y_vals, y_min)
        candidates_bottom = bottom_chain_pts[bottom_mask]
        idx_min_y = np.argmin(candidates_bottom[:, 0])  # pick smallest x if tie
        pt_min_y = candidates_bottom[idx_min_y]

        # Find point with min-x
        left_mask = np.isclose(x_vals, x_min)
        candidates_left = bottom_chain_pts[left_mask]
        idx_min_x = np.argmin(candidates_left[:, 1])  # pick smallest y if tie
        pt_min_x = candidates_left[idx_min_x]

        # Horizontal line from the min-y point to x_max
        plt.plot([pt_min_y[0], 2000], [pt_min_y[1], pt_min_y[1]],
                 color=color, linewidth=frontier_linewidth)

        # Vertical line from the min-x point to y_max
        plt.plot([pt_min_x[0], pt_min_x[0]], [pt_min_x[1], 100],
                 color=color, linewidth=frontier_linewidth)

    # ------------------------------------------------
    # 3) Plot everything
    # ------------------------------------------------
    latexify_plot(fontsize=24)
    plt.figure(figsize=(8.3, 6))

    if pareto_frontier:
        # First, scatter all points with transparency in the background
        # (No legend labels for these scatter points, so they don't clutter the legend.)
        for condition_func, color, label in filter_conditions:
            if plot_uncondensed:
                pts_unc = filtered_points[label]["uncondensed"]
                if pts_unc:
                    xs, ys = zip(*pts_unc)
                    plt.scatter(xs, ys, color=color, marker='o', alpha=point_alpha)
            if plot_condensed:
                pts_cond = filtered_points[label]["condensed"]
                if pts_cond:
                    xs, ys = zip(*pts_cond)
                    plt.scatter(xs, ys, color=color, marker='^', alpha=point_alpha)

        # Also scatter default group points
        if plot_uncondensed:
            pts_unc = default_points["uncondensed"]
            if pts_unc:
                xs, ys = zip(*pts_unc)
                plt.scatter(xs, ys, color=default_color, marker='o', alpha=point_alpha)
        if plot_condensed:
            pts_cond = default_points["condensed"]
            if pts_cond:
                xs, ys = zip(*pts_cond)
                plt.scatter(xs, ys, color=default_color, marker='^', alpha=point_alpha)

        # Now draw Pareto frontier lines
        for _, color, label in filter_conditions:
            if plot_uncondensed:
                pts_unc = filtered_points[label]["uncondensed"]
                if len(pts_unc) >= 3:
                    plot_bottom_frontier(pts_unc, color, f"{label}")
            if plot_condensed:
                pts_cond = filtered_points[label]["condensed"]
                if len(pts_cond) >= 3:
                    plot_bottom_frontier(pts_cond, color, f"{label} (condensed)")

        # Default group frontier
        if plot_uncondensed:
            pts_unc = default_points["uncondensed"]
            if len(pts_unc) >= 3:
                plot_bottom_frontier(pts_unc, default_color, f"{default_label}")
        if plot_condensed:
            pts_cond = default_points["condensed"]
            if len(pts_cond) >= 3:
                plot_bottom_frontier(pts_cond, default_color, f"{default_label} (condensed)")

    else:
        # Original scatter approach only
        for _, color, label in filter_conditions:
            if plot_uncondensed:
                pts_unc = filtered_points[label]["uncondensed"]
                if pts_unc:
                    xs, ys = zip(*pts_unc)
                    plt.scatter(xs, ys, color=color, marker='o', alpha=0.7,
                                label=f"{label}")
            if plot_condensed:
                pts_cond = filtered_points[label]["condensed"]
                if pts_cond:
                    xs, ys = zip(*pts_cond)
                    plt.scatter(xs, ys, color=color, marker='^', alpha=0.7,
                                label=f"{label} (condensed)")

        if plot_uncondensed:
            pts_unc = default_points["uncondensed"]
            if pts_unc:
                xs, ys = zip(*pts_unc)
                plt.scatter(xs, ys, color=default_color, marker='o', alpha=0.7,
                            label=f"{default_label}")
        if plot_condensed:
            pts_cond = default_points["condensed"]
            if pts_cond:
                xs, ys = zip(*pts_cond)
                plt.scatter(xs, ys, color=default_color, marker='^', alpha=0.7,
                            label=f"{default_label} (condensed)")

    # Axis labels, etc.
    if results:
        opts = results[0].get('opts', {})
        n = opts.get('n', 'n')
    else:
        n = 'n'
    
    plt.xlabel("Number of Optimization Variables")
    plt.ylabel("Increase in trace(P) [\\%]")
    plt.legend(loc='upper right', borderpad=0.2, handletextpad=0.3)
    plt.grid(True)
    #plt.tight_layout()
    plt.ylim(-1, 30)
    plt.xlim(150, 1700)
    plot_dir = get_dir("plots/linear_systems")
    filename = plot_dir / f"{filename}.pdf"
    plt.savefig(filename, bbox_inches="tight")
    plt.show()

def plot_min_lookahead_vs_optimization_vars(results, filter_conditions, trace_threshold,
                                            plot_uncondensed=True, plot_condensed=False, filename="min_line_unconstrained.pdf"):
    """
    Plots a continuous line connecting the minimal number of optimization variables
    required for a given lookahead distance for each specified class.

    Filters results based on:
    - Specified `filter_conditions`
    - A threshold for the percentage increase in trace(P)

    The line connects the **minimum** number of optimization variables found at each lookahead distance.
    """
    filtered_points = {label: {"uncondensed": {}, "condensed": {}}
                       for (_, _, label) in filter_conditions}
    
        
    # same as in plotting function above, we need to compute recompute P_LQR 
    opts = results[0].get('opts')

    m = opts["m"]
    dt = opts["dt"]
    A = results[0].get('A')
    B = results[0].get('B')

    # Hardcoded Q, R
    Q = np.eye(A.shape[0])
    R = 0.1 * np.eye(m)

    # Discretize system
    sys_gt_ct = ct.ss(A, B, np.eye(A.shape[0]), np.zeros((A.shape[0], B.shape[1])))
    sys_gt_disc_step = ct.c2d(sys_gt_ct, dt, method='zoh')
    A_d = sys_gt_disc_step.A
    B_d = sys_gt_disc_step.B

    P_LQR, _, _ = ct.dare(A_d, B_d, Q * dt, R * dt)

    for result in results:
        P_cont = result.get('P_cont')
        if P_cont is None or P_LQR is None:
            continue

        trace_P_cont = np.trace(P_cont)
        trace_P_LQR = np.trace(P_LQR)
        perc_increase = (trace_P_cont - trace_P_LQR) / trace_P_LQR * 100.0

        if abs(perc_increase) > trace_threshold or perc_increase < 0: #cancel erroneous values due to numerical isntability
            continue  # Skip results exceeding the threshold

        lookahead = result.get('t_lookahead', None)
        x_uncond = result.get('N_optim_vars_uncondensed', None)
        x_cond = result.get('N_optim_vars_condensed', None)

        if (plot_uncondensed and x_uncond is None) or (plot_condensed and x_cond is None):
            continue

        for condition_func, color, label in filter_conditions:
            if condition_func(result):
                if plot_uncondensed:
                    # Track the minimum number of uncondensed variables for each lookahead distance
                    if lookahead not in filtered_points[label]["uncondensed"]:
                        filtered_points[label]["uncondensed"][lookahead] = x_uncond
                    else:
                        filtered_points[label]["uncondensed"][lookahead] = min(
                            filtered_points[label]["uncondensed"][lookahead], x_uncond
                        )

                if plot_condensed:
                    # Track the minimum number of condensed variables for each lookahead distance
                    if lookahead not in filtered_points[label]["condensed"]:
                        filtered_points[label]["condensed"][lookahead] = x_cond
                    else:
                        filtered_points[label]["condensed"][lookahead] = min(
                            filtered_points[label]["condensed"][lookahead], x_cond
                        )

    latexify_plot(fontsize=24)
    plt.figure(figsize=(8, 6))

    # Plot the minimal number of optimization variables for each class
    for _, color, label in filter_conditions:
        if plot_uncondensed:
            sorted_lookaheads = sorted(filtered_points[label]["uncondensed"].keys())
            min_vars = [filtered_points[label]["uncondensed"][lh] for lh in sorted_lookaheads]
            if sorted_lookaheads:
                plt.plot(sorted_lookaheads, min_vars, color=color, marker='o', linestyle='-',
                         label=f"{label}")

        if plot_condensed:
            sorted_lookaheads = sorted(filtered_points[label]["condensed"].keys())
            min_vars = [filtered_points[label]["condensed"][lh] for lh in sorted_lookaheads]
            if sorted_lookaheads:
                plt.plot(sorted_lookaheads, min_vars, color=color, marker='^', linestyle='--',
                         label=f"{label} (condensed)")

    #plt.title("Minimal Number of Optimization Variables vs. Lookahead Distance")
    plt.xlabel("Time Horizon")
    plt.ylabel("Minimum Number of Variables")
    plt.legend(loc='best', fontsize=24)
    plt.grid(True)
    plt.tight_layout()
    plot_dir = get_dir("plots/linear_systems")
    filename = plot_dir / f"{filename}"
    plt.savefig(filename, bbox_inches="tight")
    plt.show()

def create_filtered_cost_vs_time_plot(results, filter_conditions,
                                      default_color='lightgray', default_label='Other Results',
                                      filter_constraints=False, epsilon=0.1,
                                      pareto_frontier=True,
                                      frontier_linewidth=2.0,
                                      point_alpha=0.3,
                                      legend=False, filename="pareto_constrained"):
    """
    Creates a plot of percentage increase in average closed-loop cost vs. mean computation time.
    
    Args:
        results (list): Loaded experiment results.
        filter_conditions (list): List of tuples (condition_function, color, label).
        default_color (str): Default color for unclassified points.
        default_label (str): Label for unclassified points.
        filter_constraints (bool): If True, only plots points that satisfy A_constr @ x <= b_constr + epsilon.
        epsilon (float): Constraint violation threshold.
        
        pareto_frontier (bool): If True, we plot a "bottom frontier" from the convex hull
                                instead of raw scatter. The raw points will still be shown
                                in the background with transparency (point_alpha).
        frontier_linewidth (float): Thickness of the frontier line and the lines from min-x and min-y.
        point_alpha (float): Transparency for points in the background when pareto_frontier=True.
    """
    # --------------------
    # 1) Organize results
    # --------------------
    filtered_points = {label: [] for (_, _, label) in filter_conditions}
    default_points = []

    # Find best (lowest) total_closed_loop_cost across all results
    all_costs = [res.get("total_closed_loop_cost", np.inf) for res in results]
    best_avg_cost = min(all_costs) if all_costs else np.inf

    for result in results:
        avg_cost = result.get("total_closed_loop_cost", None)
        mean_time = result.get("mean_solve_time", None)
        closed_loop_traj = result.get("closed_loop_trajectory", {})
        opts = result.get("experiment_opts", {})

        # Filter out NaNs
        if avg_cost is None or mean_time is None:
            continue
        if np.isnan(avg_cost) or np.isnan(mean_time):
            continue

        # Compute % increase from best
        #if best_avg_cost > 0:
        perc_increase = 100 * (avg_cost - best_avg_cost) / best_avg_cost
        # else:
        #     continue

        # Optional constraint filtering
        if filter_constraints:
            A_constr = opts.get("A_constr", None)
            b_constr = opts.get("b_constr", None)
            state_trajectory = closed_loop_traj.get("states", [])
            if A_constr is not None and b_constr is not None:
                violates_constraint = any(
                    np.any(A_constr @ x > b_constr + epsilon) for x in state_trajectory
                )
                if violates_constraint:
                    continue

        # Determine group by filter_conditions
        matched_any = False
        for condition_func, color, label in filter_conditions:
            if condition_func(result):
                filtered_points[label].append((mean_time, perc_increase, color))
                matched_any = True
                break
        if not matched_any:
            default_points.append((mean_time, perc_increase, default_color))

    # ----------------------------------
    # 2) Function to plot bottom frontier
    # ----------------------------------
    def plot_bottom_frontier(points_list, label_str):
        """
        1) Convert (x, y) into array
        2) Build convex hull
        3) Extract the bottom chain
        4) Plot it & draw lines from min-y horizontally, min-x vertically
        """
        if len(points_list) < 3:
            return

        arr = np.array([(x, y) for (x, y, c) in points_list])
        try:
            hull = ConvexHull(arr)
        except QhullError:
            return

        hull_indices = hull.vertices
        xs_hull = arr[hull_indices, 0]

        leftmost_in_hullarray = np.argmin(xs_hull)
        rightmost_in_hullarray = np.argmax(xs_hull)

        def walk_hull(start_idx, end_idx):
            chain = []
            i = start_idx
            while True:
                real_index = hull_indices[i]
                chain.append(real_index)
                if i == end_idx:
                    break
                i = (i + 1) % len(hull_indices)
            return chain

        forward_chain = walk_hull(leftmost_in_hullarray, rightmost_in_hullarray)
        backward_chain = walk_hull(rightmost_in_hullarray, leftmost_in_hullarray)

        forward_pts = arr[forward_chain, :]
        backward_pts = arr[backward_chain, :]

        if np.mean(forward_pts[:, 1]) < np.mean(backward_pts[:, 1]):
            bottom_chain_pts = forward_pts
        else:
            bottom_chain_pts = backward_pts

        color_use = points_list[0][2]  # color from the first triple (x,y,color)
        # If label_str is None, don't pass a label to avoid a second legend entry
        plt.plot(bottom_chain_pts[:, 0], bottom_chain_pts[:, 1],
                 color=color_use, linewidth=frontier_linewidth,
                 label=label_str if label_str is not None else None)

        x_vals = bottom_chain_pts[:, 0]
        y_vals = bottom_chain_pts[:, 1]
        x_min = np.min(x_vals)
        x_max = np.max(x_vals)
        y_min = np.min(y_vals)
        y_max = np.max(y_vals)

        bottom_mask = np.isclose(y_vals, y_min)
        bottom_candidates = bottom_chain_pts[bottom_mask]
        idx_min_y = np.argmin(bottom_candidates[:, 0])
        pt_min_y = bottom_candidates[idx_min_y]

        left_mask = np.isclose(x_vals, x_min)
        left_candidates = bottom_chain_pts[left_mask]
        idx_min_x = np.argmin(left_candidates[:, 1])
        pt_min_x = left_candidates[idx_min_x]

        # These lines have fixed endpoints (user-specified)
        plt.plot([pt_min_y[0], 0.01],
                 [pt_min_y[1], pt_min_y[1]],
                 color=color_use, linewidth=frontier_linewidth)
        plt.plot([pt_min_x[0], pt_min_x[0]],
                 [pt_min_x[1], 100],
                 color=color_use, linewidth=frontier_linewidth)

    # ----------------------------------
    # 3) Plot
    # ----------------------------------
    latexify_plot(fontsize=28)
    plt.figure(figsize=(9.4, 6))

    if pareto_frontier:
        # Scatter points in the background with alpha=point_alpha
        for label, pts in filtered_points.items():
            if pts:
                xs = [p[0] for p in pts]
                ys = [p[1] for p in pts]
                colors = [p[2] for p in pts]
                c0 = colors[0]

                # 'Fixed Model (condensed)' => bigger crosses + legend, in foreground
                if label == 'Fixed Model (condensed)':
                    plt.scatter(xs, ys, color=c0, marker='x', s=200,
                                alpha=1.0, label=label, zorder=3)
                else:
                    # For other classes, just dots, no legend
                    plt.scatter(xs, ys, color=c0, marker='o', alpha=point_alpha)

        if default_points:
            xs = [p[0] for p in default_points]
            ys = [p[1] for p in default_points]
            plt.scatter(xs, ys, color=default_color, marker='o', alpha=point_alpha)

        # Now plot the Pareto frontier lines
        for (condition_func, color, label) in filter_conditions:
            pts = filtered_points[label]
            if len(pts) >= 3:
                # For 'Fixed Model (condensed)' => no second legend entry
                if label == 'Fixed Model (condensed)':
                    plot_bottom_frontier(pts, None)
                else:
                    plot_bottom_frontier(pts, label)

        if len(default_points) >= 3:
            plot_bottom_frontier(default_points, default_label)

    else:
        # Original approach: scatter + legend
        for label, points in filtered_points.items():
            if points:
                xs, ys, colors = zip(*points)
                c0 = colors[0]

                if label == 'Fixed Model (condensed)':
                    plt.scatter(xs, ys, color=c0, marker='x', s=100,
                                alpha=0.7, label=label, zorder=3)
                else:
                    plt.scatter(xs, ys, color=c0, marker='o', alpha=0.7, label=label)

        if default_points:
            xs, ys, _ = zip(*default_points)
            plt.scatter(xs, ys, color=default_color, marker='o', alpha=0.7, label=default_label)

    if results:
        opts = results[0].get('experiment_opts', {})
        n = opts.get('n', 'n')
    else:
        n = 'n'

    #plt.title(f"% Increase in Avg. Closed-Loop Cost vs. Mean Solve Time (n={n})")
    plt.xlabel("Mean Solve Time [s]")
    plt.ylabel("Increase in Closed-Loop Cost [\%]")
    plt.ylim(-1, 30)
    plt.xlim(0, 0.006)
    if legend:
        plt.legend(loc='best', fontsize=28)
    plt.grid(True)
    plt.tight_layout()
    plot_dir = get_dir("plots/linear_systems")
    filename = plot_dir / f"{filename}.pdf"
    plt.savefig(filename, bbox_inches="tight")
    plt.show()