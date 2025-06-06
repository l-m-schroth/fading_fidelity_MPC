import numpy as np
import matplotlib.pyplot as plt
import mujoco
import os
import mediapy as media
import casadi as ca
from copy import deepcopy
from scipy.optimize import fsolve

def generate_exponential_step_sizes(initial_step: float, rate: float, num_steps: int, plot: bool = False) -> list:
    """
    Generate a list of exponentially increasing step sizes.

    Parameters:
        initial_step (float): The initial step size.
        rate (float): The exponential growth rate (e.g., 1.05 for a 5% increase per step).
        num_steps (int): The total number of step sizes to generate.
        plot (bool): If True, display a plot of the step size development.

    Returns:
        list: A list of exponentially increasing step sizes.
    """
    step_sizes = [initial_step * (rate ** i) for i in range(num_steps)]

    if plot:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(range(num_steps), step_sizes, marker='o')
        plt.title("Exponential Step Sizes")
        plt.xlabel("Step Index")
        plt.ylabel("Step Size")
        plt.grid(True)
        plt.show()

    return step_sizes

def get_x_and_y_pos(angles, link_length):
    # Initialize the starting point (top of the chain)
    x, y = [0], [0]  # Starting at the origin (0, 0)
    
    # Initialize cumulative angle, starting downward
    cumulative_angle = 0 # Downward direction
    
    for angle in angles:
        # Update cumulative angle (relative to the downward direction)
        cumulative_angle += angle
        
        # Calculate next joint position
        x_next = x[-1] + link_length * np.sin(cumulative_angle)
        y_next = y[-1] + -link_length * np.cos(cumulative_angle)
        
        # Append the new joint position
        x.append(x_next)
        y.append(y_next)

    return x, y  


def get_ee_position(data, ee_site_id, steady_state_z_values):
    return data.site_xpos[ee_site_id] - np.array([0, 0, steady_state_z_values[0]])

# Function to reset simulation after initialization
def reset_sim(mjData, qpos, qvel):
    mjData.qpos[:] = qpos
    mjData.qvel[:] = qvel
    mjData.time = 0.0

# Function to check convergence to the origin
def has_converged(states, threshold, window_size):
    if len(states) < window_size:
        return False
    window = np.array(states[-window_size:])
    return np.all(np.abs(window) < threshold)

def get_ee_position(data, ee_site_id, steady_state_z_values):
    return data.site_xpos[ee_site_id] - np.array([0, 0, steady_state_z_values[0]])  

def get_x_and_y_pos_casadi(angles, link_length):
    # Initialize symbolic variables
    x = [0]
    y = [0]
    cumulative_angle = 0  # Starting downward

    # Use CasADi's vertsplit to iterate over symbolic elements
    for angle in ca.vertsplit(angles):
        cumulative_angle += angle
        x_next = x[-1] + link_length * ca.sin(cumulative_angle)
        y_next = y[-1] + -link_length * ca.cos(cumulative_angle)
        x.append(x_next)
        y.append(y_next)

    return ca.vertcat(*x), ca.vertcat(*y)

def plot_chains(link_lengths, angles_list, fontsizes, save=False):
    plt.figure(figsize=(8, 8))
    if not save:
        plt.title("Hanging 2D Rigid Link Chains", fontsize=fontsizes["title"])
    plt.xlabel("X-axis", fontsize=fontsizes["labels"])
    plt.ylabel("Y-axis", fontsize=fontsizes["labels"])
    plt.grid(True)
    plt.axis('equal')  # Equal aspect ratio for correct visualization

    for link_length, angles in zip(link_lengths, angles_list):
        x, y = get_x_and_y_pos(angles, link_length)
        
        # Plot the chain
        num_links = len(angles)
        plt.plot(x, y, marker='o', linestyle='-', linewidth=2, markersize=8, label=f"Chain with {num_links} links")

    plt.legend(fontsize=fontsizes["legend"])
    plt.xticks(fontsize=fontsizes["ticks"])
    plt.yticks(fontsize=fontsizes["ticks"])
    if save:
        os.makedirs("plots", exist_ok=True)
        filename = f"plots/chain_configurations.pgf"
        plt.savefig(filename, bbox_inches="tight")
    plt.show()

def process_list(link_length, x_prev_list, y_prev_list):
    x_list, y_list = [0], [0]
    angles = []
    prev_angle_tot = 0
    num_links = int((len(x_prev_list)-1)/2)
    for i in range(1,num_links+1):
        x_prev = x_prev_list[2*i]
        y_prev = y_prev_list[2*i]
        x_diff =  x_prev - x_list[-1]
        y_diff =  y_list[-1] - y_prev
        next_angle_tot = np.arctan2(x_diff, y_diff) # for AD later: atan2 is discontinuous, but the derivative can be continously defined except for the origin -> x_diff=0 and y_diff=0 does not happen, so we should be good.
        x_next = x_list[-1] + np.sin(next_angle_tot)*link_length
        y_next = y_list[-1] - np.cos(next_angle_tot)*link_length
        
        if i > 1:
            next_angle_rel = next_angle_tot - prev_angle_tot
        else:
            next_angle_rel = next_angle_tot

        # update previous total angle
        prev_angle_tot = next_angle_tot

        # append to lists
        angles.append(next_angle_rel)
        x_list.append(x_next)
        y_list.append(y_next)

    return angles, x_list, y_list

# --- fucntion to heuristically map the joint angles to joint angles of lower dimensional chains --- #
def compute_q(q_high_dim, link_length_high_dim, link_length_low_dim, n_chains_low):
    n_chains_high = len(q_high_dim)
    x_high, y_high = get_x_and_y_pos(q_high_dim, link_length_high_dim)
    step_size = int(n_chains_high/n_chains_low)
    q_low = []
    q_tot, x_actual, y_actual = 0, 0, 0
    for i in range(n_chains_low):
        x_diff = x_high[step_size*(i+1)] - x_actual
        y_diff = y_high[step_size*(i+1)] - y_actual
        q_low_next = np.arctan2(x_diff, -y_diff) - q_tot
        q_low.append(q_low_next)
        q_tot += q_low_next
        x_actual += np.sin(q_tot)*link_length_low_dim
        y_actual += -np.cos(q_tot)*link_length_low_dim
    return q_low

def compute_q_casadi(q_high_dim, link_length_high_dim, link_length_low_dim, n_chains_low):
    # Get dimensions
    n_chains_high = q_high_dim.size()[0]  # Size in CasADi syntax
    step_size = n_chains_high // n_chains_low

    # Symbolic variables
    x_high, y_high = get_x_and_y_pos_casadi(q_high_dim, link_length_high_dim)
    
    # Initialize low-dimensional chain computation
    q_low = []
    q_tot = 0
    x_actual, y_actual = 0, 0
    
    for i in range(n_chains_low):
        x_diff = x_high[step_size * (i + 1)] - x_actual
        y_diff = y_high[step_size * (i + 1)] - y_actual
        q_low_next = ca.arctan2(x_diff, -y_diff) - q_tot
        q_low.append(q_low_next)
        q_tot += q_low_next
        x_actual += ca.sin(q_tot) * link_length_low_dim
        y_actual += -ca.cos(q_tot) * link_length_low_dim

    return ca.vertcat(*q_low)  # Return as CasADi symbolic vector

# unified compute controls function that can handle all cases
def compute_controls(u_1_high_dim, u_2_high_dim, link_length_high_dim, link_length_low_dim, stiffness_high_dim, stiffness_low_dim, gear, n_chains_high, n_chains_low):
    n_links_per_angle_high_dim = n_chains_high // 2
    n_links_per_angle_low_dim = n_chains_low // 2
    q_high_1 = gear*u_1_high_dim/stiffness_high_dim
    q_high_2 = gear*u_2_high_dim/stiffness_high_dim
    x_high, y_high = get_x_and_y_pos([q_high_1]*n_links_per_angle_high_dim + [q_high_2]*n_links_per_angle_high_dim, link_length_high_dim)
    q_1_desired_jointly = np.arctan2(x_high[n_links_per_angle_high_dim], -y_high[n_links_per_angle_high_dim])
    angle_factor = 2 / (n_links_per_angle_low_dim + 1)
    q_1 = angle_factor * q_1_desired_jointly
    x_actual = 0
    y_actual = 0
    for i in range(1, n_links_per_angle_low_dim + 1):
        x_actual += np.sin(i*q_1)*link_length_low_dim
        y_actual -= np.cos(i*q_1)*link_length_low_dim
    x_diff = x_high[-1] - x_actual
    y_diff = y_high[-1] - y_actual
    q_2_desired_jointly = np.arctan2(x_diff, -y_diff) - n_links_per_angle_low_dim*q_1
    q_2 = angle_factor*q_2_desired_jointly
    u_1 = q_1 * stiffness_low_dim / gear
    u_2 = q_2 * stiffness_low_dim / gear   
    return np.array([u_1, u_2])


def compute_controls_inverse(u_1_low_dim, u_2_low_dim, link_length_high_dim, link_length_low_dim, stiffness_high_dim, stiffness_low_dim, gear, n_chains_high, n_chains_low, heuristic=True):
    n_links_per_angle_high_dim = n_chains_high // 2
    n_links_per_angle_low_dim = n_chains_low // 2

    q_low_1 = gear * u_1_low_dim / stiffness_low_dim
    q_low_2 = gear * u_2_low_dim / stiffness_low_dim
    
    # Assuming get_x_and_y_pos is defined elsewhere
    x_low_dim, y_low_dim = get_x_and_y_pos(
        [q_low_1] * n_links_per_angle_low_dim + [q_low_2] * n_links_per_angle_low_dim,
        link_length_low_dim
    )
    q_1_desired_jointly = np.arctan2(x_low_dim[n_links_per_angle_low_dim], -y_low_dim[n_links_per_angle_low_dim])

    if not heuristic:
        # Define the equation for root finding for q_1
        def q1_equation(q):
            x_high_dim = 0
            y_high_dim = 0
            for i in range(1, n_links_per_angle_high_dim + 1):
                x_high_dim += np.sin(i * q)
                y_high_dim -= np.cos(i * q)
            return np.tan(q_1_desired_jointly) - x_high_dim / y_high_dim

        # Solve for q_1
        q_1 = fsolve(q1_equation, q_1_desired_jointly * 2 / (n_links_per_angle_high_dim + 1))[0]
    else:
        q_1 = q_1_desired_jointly * 2 / (n_links_per_angle_high_dim + 1)

    x_high_dim_actual = 0
    y_high_dim_actual = 0
    for i in range(1, n_links_per_angle_high_dim + 1):
        x_high_dim_actual += np.sin(i * q_1)*link_length_high_dim
        y_high_dim_actual -= np.cos(i * q_1)*link_length_high_dim

    x_diff = x_low_dim[-1] - x_high_dim_actual
    y_diff = y_low_dim[-1] - y_high_dim_actual
    q_2_desired_jointly = np.arctan2(x_diff, -y_diff) - n_links_per_angle_high_dim * q_1

    if not heuristic:
        # Define the equation for root finding for q_2
        def q2_equation(q):
            x_high_dim = 0
            y_high_dim = 0
            for i in range(1, n_links_per_angle_high_dim + 1):
                x_high_dim += np.sin(i * q)
                y_high_dim -= np.cos(i * q)
            return np.tan(q_2_desired_jointly) - x_high_dim / y_high_dim

        # Solve for q_2
        q_2 = fsolve(q2_equation, q_2_desired_jointly * 2 / (n_links_per_angle_high_dim + 1))[0]
    else:
        q_2 = q_2_desired_jointly * 2 / (n_links_per_angle_high_dim + 1)

    u_1 = q_1 * stiffness_high_dim / gear
    u_2 = q_2 * stiffness_high_dim / gear
    
    return np.array([u_1, u_2])


def compute_controls_casadi(u_1_high_dim, u_2_high_dim, link_length_high_dim, link_length_low_dim, stiffness_high_dim, stiffness_low_dim, gear, n_chains_high, n_chains_low):
    n_links_per_angle_high_dim = n_chains_high // 2
    n_links_per_angle_low_dim = n_chains_low // 2
    q_high_1 = gear * u_1_high_dim / stiffness_high_dim
    q_high_2 = gear * u_2_high_dim / stiffness_high_dim
    x_high, y_high = get_x_and_y_pos_casadi(ca.vertcat(*([q_high_1] * n_links_per_angle_high_dim + [q_high_2] * n_links_per_angle_high_dim)), link_length_high_dim)
    q_1_desired_jointly = ca.arctan2(x_high[n_links_per_angle_high_dim], -y_high[n_links_per_angle_high_dim])
    angle_factor = 2 / (n_links_per_angle_low_dim + 1)
    q_1 = angle_factor * q_1_desired_jointly
    x_actual = 0
    y_actual = 0
    for i in range(1, n_links_per_angle_low_dim + 1):
        x_actual += ca.sin(i * q_1) * link_length_low_dim
        y_actual -= ca.cos(i * q_1) * link_length_low_dim
    x_diff = x_high[-1] - x_actual
    y_diff = y_high[-1] - y_actual
    q_2_desired_jointly = ca.arctan2(x_diff, -y_diff) - n_links_per_angle_low_dim * q_1
    q_2 = angle_factor * q_2_desired_jointly
    u_1 = q_1 * stiffness_low_dim / gear
    u_2 = q_2 * stiffness_low_dim / gear
    return ca.vertcat(u_1, u_2)

def clip_angle(angle):
    """Clip angles to the range (-π, π]."""
    return (angle + np.pi) % (2 * np.pi) - np.pi

def compute_controls_smooth(u_1_high_dim, u_2_high_dim, link_length_high_dim, link_length_low_dim, stiffness_high_dim, stiffness_low_dim, gear, n_chains_high, n_chains_low):
    n_links_per_angle_high_dim = n_chains_high // 2
    n_links_per_angle_low_dim = n_chains_low // 2
    q_high_1 = gear * u_1_high_dim / stiffness_high_dim
    q_high_2 = gear * u_2_high_dim / stiffness_high_dim
    
    x_high, y_high = get_x_and_y_pos(
        [clip_angle(q_high_1)] * n_links_per_angle_high_dim + [clip_angle(q_high_2)] * n_links_per_angle_high_dim,
        link_length_high_dim
    )
    
    q_1_desired_jointly = clip_angle(np.arctan2(x_high[n_links_per_angle_high_dim], -y_high[n_links_per_angle_high_dim]))
    angle_factor = 2 / (n_links_per_angle_low_dim + 1)
    q_1 = angle_factor * clip_angle(q_1_desired_jointly)
    
    x_actual = 0
    y_actual = 0
    for i in range(1, n_links_per_angle_low_dim + 1):
        x_actual += np.sin(i * q_1) * link_length_low_dim
        y_actual -= np.cos(i * q_1) * link_length_low_dim
    
    x_diff = x_high[-1] - x_actual
    y_diff = y_high[-1] - y_actual
    q_2_desired_jointly = clip_angle(np.arctan2(x_diff, -y_diff) - n_links_per_angle_low_dim * q_1)
    q_2 = angle_factor * q_2_desired_jointly
    
    u_1 = q_1 * stiffness_low_dim / gear
    u_2 = q_2 * stiffness_low_dim / gear   
    return np.array([u_1, u_2])

def compute_controls_old(u_1_16, u_2_16, link_length_list, stiffness_list, gear):
    # old compute controls function, not used anymore, better version
    u_controls = [np.array([u_1_16, u_2_16])]
    stiffness_16 = stiffness_list[0]
    q_16_1 = gear*u_1_16/stiffness_16
    q_16_2 = gear*u_2_16/stiffness_16
    q_16 = [q_16_1]*8 + [q_16_2]*8
    x_prev, y_prev = get_x_and_y_pos(q_16, link_length_list[0]) # get positions from the 16 link chain
    for i in range(len(link_length_list)-1):
        q_next, x_next, y_next = process_list(link_length_list[i+1], x_prev, y_prev)
        middle = int(len(q_next)/2)
        match i:
            case 0:
                # 8 links case
                q_1 = (4*q_next[0] + 3*q_next[1] + 2*q_next[2] + q_next[3])/10
                q_2 = (4*q_next[4] + 3*q_next[5] + 2*q_next[6] + q_next[7])/10
            case 1:
                # 4 link case
                q_1 = (2*q_next[0]+q_next[1])/3
                q_2 = (2*q_next[2]+q_next[3])/3
            case 2:
                # 2 link case  
                q_1 = q_next[0] 
                q_2 = q_next[1] 
        u_1 = q_1 * stiffness_list[i+1] / gear
        u_2 = q_2 * stiffness_list[i+1] / gear        
        u_controls.append(np.array([u_1, u_2]))
        # assign new values
        x_prev, y_prev = x_next, y_next
    return u_controls

def simulate(mjModel, mjData, duration, ee_site_id, qpos_init, qvel_init, create_video=False, u_controls = None):

    # reset sim to inital condition
    reset_sim(mjData, qpos_init, qvel_init)

    # init lists
    mujoco.mj_fwdPosition(mjModel, mjData)
    ee_positions = []
    q_positions = [qpos_init] 
    q_velocities = [qvel_init] 
    u_counter = 0
    if u_controls is not None:
        if isinstance(u_controls, list):
            u_controls = np.array(u_controls)
        if len(u_controls.shape) > 1:
            multiple_u = True
        else:
            multiple_u = False
    if create_video:
        frames = []
        framerate = 60

    with mujoco.Renderer(mjModel) as renderer:
        while mjData.time <= duration - mjModel.opt.timestep:

            if u_controls is not None:
                middle = int(len(mjData.ctrl[:])/2)
                if multiple_u: 
                    u = u_controls[u_counter, :]
                else:
                    u = u_controls    
                mjData.ctrl[:middle] = u[0]
                mjData.ctrl[middle:] = u[1]
                u_counter += 1

            mujoco.mj_step(mjModel, mjData)
            
            ee_pos = get_ee_position(mjData, ee_site_id, [0.2]) # ee steady state z_value is 0.2
            ee_positions.append(ee_pos)

            # Render the frame
            if create_video and len(frames) < mjData.time * framerate:
                renderer.update_scene(mjData)
                pixels = renderer.render()
                frames.append(pixels)

            #print(mjData.qpos)
            #print(mjData.qvel)
            q_positions.append(deepcopy(mjData.qpos))
            q_velocities.append(deepcopy(mjData.qvel))

    # Save rendered video
    if create_video:
        current_dir = os.path.dirname(__file__)
        video_folder = os.path.join(current_dir, "trajectories", "videos")
        if u_controls is not None:
            output_video_file = os.path.join(video_folder, f"video_chain_comparison_{len(qpos_init)}_controlled.mp4")
        else:    
            output_video_file = os.path.join(video_folder, f"video_chain_comparison_{len(qpos_init)}.mp4")
        media.write_video(output_video_file, frames, fps=framerate)
        print(f"Simulation video saved to {output_video_file}.")

    return ee_positions, q_positions, q_velocities