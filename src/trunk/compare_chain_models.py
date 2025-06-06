"""
This file provides additional visualizations for the chain models and the mapping between the models. Also things like the autonomous behaviour is visualized.
"""
import mujoco
import numpy as np
import os
import matplotlib.pyplot as plt
import casadi as ca
from matplotlib.animation import FuncAnimation
from copy import deepcopy
from acados_template import latexify_plot

from trunk.trunk_utils import get_x_and_y_pos, process_list, plot_chains, simulate, compute_controls, compute_q_casadi, compute_q

# NOTE: Currently the starting angles are computed based on the positions of the PREVIOUS chain instead of the chain with the most links, which would maybe be a bit more accurate, but it already looks good visually. 
# Could investigate if it makes a difference

print(mujoco.__version__)

current_dir = os.path.dirname(__file__)

create_video = False
map_controls = False # maps controls, currently switching off only affects sinusoidal control

# fontsizes for plots, some plots still need adaption, not all make use of the struct
latexify_plot(fontsize=18)
fontsizes = {}
fontsizes["labels"] = 18#15
fontsizes["legend"] = 18#15
fontsizes["title"] = 18#17
fontsizes["ticks"] = 18#15

# load 2 link model
xml_path_2 = os.path.join(current_dir, "models", "chain_models", "chain_2_links_expanded.xml")
model_2 = mujoco.MjModel.from_xml_path(xml_path_2)
data_2 = mujoco.MjData(model_2)
qpos_2_ss = data_2.qpos

# load 4 link model
xml_path_4 = os.path.join(current_dir, "models", "chain_models", "chain_4_links_expanded.xml")
model_4 = mujoco.MjModel.from_xml_path(xml_path_4)
data_4 = mujoco.MjData(model_4)
qpos_4_ss = data_4.qpos

# load 8 link model
xml_path_8 = os.path.join(current_dir, "models", "chain_models", "chain_8_links_expanded.xml")
model_8 = mujoco.MjModel.from_xml_path(xml_path_8)
data_8 = mujoco.MjData(model_8)
qpos_8_ss = data_8.qpos

# load 16 link model
xml_path_16 = os.path.join(current_dir, "models", "chain_models", "chain_16_links_expanded.xml")
model_16 = mujoco.MjModel.from_xml_path(xml_path_16)
data_16 = mujoco.MjData(model_16)
qpos_16_ss = data_16.qpos

# get site IDs for all models
ee_site_id_2 = mujoco.mj_name2id(model_2, mujoco.mjtObj.mjOBJ_SITE, "endeffector")
ee_site_id_4 = mujoco.mj_name2id(model_4, mujoco.mjtObj.mjOBJ_SITE, "endeffector")
ee_site_id_8 = mujoco.mj_name2id(model_8, mujoco.mjtObj.mjOBJ_SITE, "endeffector")
ee_site_id_16 = mujoco.mj_name2id(model_16, mujoco.mjtObj.mjOBJ_SITE, "endeffector")

# --- initial angles and angular velocities ---
angle_1 = np.pi/20
angle_2 = -np.pi/45

# 16 link chain
link_length_16 = 0.03125
qpos_16_init = [angle_1]*8 + [angle_2]*8 
# let's start by setting the velocity to zero 
qvel_16_init = [0]*16
x_16, y_16 = get_x_and_y_pos(qpos_16_init, link_length_16)

# 8 link chain
link_length_8 = link_length_16*2
qpos_8_init, x_8, y_8 = process_list(link_length_8, x_16, y_16)
qvel_8_init = [0]*8

# 4 link chain
link_length_4 = link_length_8*2
qpos_4_init, x_4, y_4 = process_list(link_length_4, x_8, y_8)
qvel_4_init = [0]*4

# 2 link chain
link_length_2 = link_length_4*2
qpos_2_init, x_2, y_2 = process_list(link_length_2, x_4, y_4)
qvel_2_init = [0]*2

#verify that everything is as intended
link_length_list = [link_length_16, link_length_8, link_length_4, link_length_2]
qpos_list = [qpos_16_init, qpos_8_init, qpos_4_init, qpos_2_init]
plot_chains(link_length_list, qpos_list, fontsizes, save=True)

# simulate systems
duration = 6.0

ee_positions_16, q_positions_16, _ = simulate(model_16, deepcopy(data_16), duration, ee_site_id_16, qpos_16_init, qvel_16_init, create_video)
ee_positions_8, q_positions_8, _ = simulate(model_8, deepcopy(data_8), duration, ee_site_id_8, qpos_8_init, qvel_8_init, create_video)
ee_positions_4, q_positions_4, _ = simulate(model_4, deepcopy(data_4), duration, ee_site_id_4, qpos_4_init, qvel_4_init, create_video)
ee_positions_2, q_positions_2, _ = simulate(model_2, deepcopy(data_2), duration, ee_site_id_2, qpos_2_init, qvel_2_init, create_video)

# plotting

# Extract x and z positions over time for each system
time = np.linspace(0, duration, len(ee_positions_16))

ee_positions_16_x = [pos[0] for pos in ee_positions_16]
ee_positions_16_z = [pos[2] for pos in ee_positions_16]

ee_positions_8_x = [pos[0] for pos in ee_positions_8]
ee_positions_8_z = [pos[2] for pos in ee_positions_8]

ee_positions_4_x = [pos[0] for pos in ee_positions_4]
ee_positions_4_z = [pos[2] for pos in ee_positions_4]

ee_positions_2_x = [pos[0] for pos in ee_positions_2]
ee_positions_2_z = [pos[2] for pos in ee_positions_2]

# Plot endeffector x and z positions over time in a 2x1 subplot
end_index = 500
fig, axes = plt.subplots(2, 1, figsize=(10, 8))
axes[0].plot(time[:end_index], ee_positions_16_x[:end_index], label="16-link")
axes[0].plot(time[:end_index], ee_positions_8_x[:end_index], label="8-link")
axes[0].plot(time[:end_index], ee_positions_4_x[:end_index], label="4-link")
axes[0].plot(time[:end_index], ee_positions_2_x[:end_index], label="2-link")
#axes[0].set_title('Endeffector X Position Over Time')
axes[0].set_xlabel('Time (s)')
axes[0].set_ylabel('x')
axes[0].grid(True)
axes[0].legend()

axes[1].plot(time[:end_index], ee_positions_16_z[:end_index], label="16-link")
axes[1].plot(time[:end_index], ee_positions_8_z[:end_index], label="8-link")
axes[1].plot(time[:end_index], ee_positions_4_z[:end_index], label="4-link")
axes[1].plot(time[:end_index], ee_positions_2_z[:end_index], label="2-link")
#axes[1].set_title('Endeffector Z Position Over Time')
axes[1].set_xlabel('Time (s)')
axes[1].set_ylabel('z')
axes[1].grid(True)

#axes[1].legend()

plt.tight_layout()
plt.show()

# Plot development in the x-z plane
plt.figure(figsize=(10, 6))
plt.plot(ee_positions_16_x, ee_positions_16_z, label="16-link")
plt.plot(ee_positions_8_x, ee_positions_8_z, label="8-link")
plt.plot(ee_positions_4_x, ee_positions_4_z, label="4-link")
plt.plot(ee_positions_2_x, ee_positions_2_z, label="2-link")
#plt.title('Endeffector Position in X-Z Plane')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid()
plt.show()

def animate_chains(q_trajectories, link_lengths, duration, fontsizes, interval=50):
    """
    Creates an animated plot showing the progression of multiple chains over time.

    Parameters:
    - q_trajectories: List of numpy arrays containing joint trajectories for each chain.
    - link_lengths: List of link lengths for each chain.
    - duration: Duration of the animation in seconds.
    - interval: Interval between frames in milliseconds (default is 50ms).
    """
    num_chains = len(q_trajectories)
    colors = ["blue", "green", "red", "purple"]  # Colors for the chains
    time_steps = len(q_trajectories[0])

    # Initialize the figure
    fig, ax = plt.subplots(figsize=(8, 6))
    lines = []
    for i in range(num_chains):
        line, = ax.plot([], [], lw=2, color=colors[i % len(colors)], label=f"{2**(i+1)}-link chain")
        lines.append(line)

    ax.set_xlim(-0.2, 0.4)  # Adjust as needed for the visualization
    ax.set_ylim(-0.6, 0.1)
    ax.set_xlabel("X Position", fontsize=fontsizes["labels"])
    ax.set_ylabel("Y Position", fontsize=fontsizes["labels"])
    ax.set_title("Autonomous Behaviour of the Chains", fontsize=fontsizes["title"])
    ax.legend(fontsize=fontsizes["legend"])
    ax.tick_params(axis='x', labelsize=fontsizes["ticks"])  # x-axis ticks
    ax.tick_params(axis='y', labelsize=fontsizes["ticks"])  # y-axis ticks
    ax.grid()
    
    # Compute time values
    time = np.linspace(0, duration, time_steps)

    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    def update(frame):
        for i, q in enumerate(q_trajectories):
            x, y = get_x_and_y_pos(q[frame], link_lengths[i])
            lines[i].set_data(x, y)
        return lines

    ani = FuncAnimation(
        fig, update, frames=time_steps, init_func=init, blit=True, interval=interval
    )

    plt.show()

    return ani


# Call the function to animate the chains
ani = animate_chains(
    q_trajectories=[q_positions_16, q_positions_8, q_positions_4, q_positions_2],
    link_lengths=[link_length_16, link_length_8, link_length_4, link_length_2],
    fontsizes=fontsizes,
    duration=duration
)


# --- compare behaviour for actuated trajectories ---

u_1_16 = 0.1 # controls in the first section
u_2_16 = -0.5 # controls in the second section

stiffness_list = [3, 1.55, 0.9, 0.65]
stiffness_16 = 3
stiffness_8 = 1.55
stiffness_4 = 0.9
stiffness_2 = 0.65
gear = 2

# u_controls_16, u_controls_8, u_controls_4, u_controls_2 = compute_controls(u_1, u_2, link_length_list, stiffness_list, gear)
u_controls_2 = compute_controls(u_1_16, u_2_16, link_length_16, link_length_2, stiffness_16, stiffness_2, gear, 16, 2)
u_controls_4 = compute_controls(u_1_16, u_2_16, link_length_16, link_length_4, stiffness_16, stiffness_4, gear, 16, 4)
u_controls_8 = compute_controls(u_1_16, u_2_16, link_length_16, link_length_8, stiffness_16, stiffness_8, gear, 16, 8)
u_controls_16 = np.array([u_1_16, u_2_16])

duration_controlled = 3.0

if not map_controls:
    # try overwriting just to see what happens REMOVE LATER
    u_controls_2 = u_controls_16
    u_controls_4 = u_controls_16
    u_controls_8 = u_controls_16

data_16_inp = deepcopy(data_16)
data_8_inp = deepcopy(data_8)
data_4_inp = deepcopy(data_4)
data_2_inp = deepcopy(data_2)

# simulate with controls
ee_positions_16_controlled, qpos_16_controlled, _ = simulate(model_16, data_16_inp, duration_controlled, ee_site_id_16, qpos_16_ss, qvel_16_init, create_video, u_controls_16)
ee_positions_8_controlled, qpos_8_controlled, _ = simulate(model_8, data_8_inp, duration_controlled, ee_site_id_8, qpos_8_ss, qvel_8_init, create_video, u_controls_8)
ee_positions_4_controlled, qpos_4_controlled, _ = simulate(model_4, data_4_inp, duration_controlled, ee_site_id_4, qpos_4_ss, qvel_4_init, create_video, u_controls_4)
ee_positions_2_controlled, qpos_2_controlled, _ = simulate(model_2, data_2_inp, duration_controlled, ee_site_id_2, qpos_2_ss, qvel_2_init, create_video, u_controls_2)

# Obtain final resulting configurations 
final_qpos_16 = qpos_16_controlled[-1]
final_qpos_8 = qpos_8_controlled[-1]
final_qpos_4 = qpos_4_controlled[-1]
final_qpos_2 = qpos_2_controlled[-1]
final_qvel_16 = data_16.qvel
final_qvel_8 = data_8.qvel
final_qvel_4 = data_4.qvel
final_qvel_2 = data_2.qvel

plot_chains([link_length_16, link_length_8, link_length_4, link_length_2], [final_qpos_16, final_qpos_8, final_qpos_4, final_qpos_2], fontsizes)


# # Compute the actually desired first angle
# x_16_controlled, y_16_controlled = get_x_and_y_pos(final_qpos_16, link_length_16)
# qpos_8_controlled_desired, x_8_controlled, y_8_controlled = process_list(link_length_8, x_16_controlled, y_16_controlled)
# qpos_4_controlled_desired, x_4_controlled, y_4_controlled = process_list(link_length_4, x_8_controlled, y_8_controlled)
# qpos_2_controlled_desired, x_2_controlled, y_2_controlled = process_list(link_length_2, x_4_controlled, y_4_controlled)

# # print ratios between desired and actual angle 
# ratios_8 = qpos_8_controlled_desired/final_qpos_8
# ratios_4 = qpos_4_controlled_desired/final_qpos_4
# ratios_2 =  qpos_2_controlled_desired/final_qpos_2

# # print ratios for tuning, scale the gear ratio of the actautors correspondingly, want to achieve 1 here 
# print("Ratios for 8 link chain:", ratios_8)
# print("Mean ratio for 8 link chain:", np.mean(ratios_8))
# print("Ratios for 4 link chain:", ratios_4)
# print("Mean ratios for 4 link chain:", np.mean(ratios_4))
# print("Ratios for 2 link chain:", ratios_2)
# print("Mean ratios for 2 link chain:", np.mean(ratios_2))

# --- try switching the control input after some time ---
switch_input = True
if switch_input:
    u_1_16_2 = 0.3 # controls in the first section
    u_2_16_2 = 0.3 # controls in the second section
    create_video = False

    u_controls_2_2 = compute_controls(u_1_16_2, u_2_16_2, link_length_16, link_length_2, stiffness_16, stiffness_2, gear, 16, 2)
    u_controls_4_2 = compute_controls(u_1_16_2, u_2_16_2, link_length_16, link_length_4, stiffness_16, stiffness_4, gear, 16, 4)
    u_controls_8_2 = compute_controls(u_1_16_2, u_2_16_2, link_length_16, link_length_8, stiffness_16, stiffness_8, gear, 16, 8)
    u_controls_16_2 = np.array([u_1_16_2, u_2_16_2])

    duration_controlled_2 = 4.0

    if not map_controls:
        # try overwriting just to see what happens REMOVE LATER
        u_controls_2_2 = u_controls_16_2
        u_controls_4_2 = u_controls_16_2
        u_controls_8_2 = u_controls_16_2

    # simulate with controls
    ee_positions_16_controlled_2, qpos_16_controlled_2, _ = simulate(model_16, data_16_inp, duration_controlled_2, ee_site_id_16, final_qpos_16, final_qvel_16, create_video, u_controls_16_2)
    ee_positions_8_controlled_2, qpos_8_controlled_2, _ = simulate(model_8, data_8_inp, duration_controlled_2, ee_site_id_8, final_qpos_8, final_qvel_8, create_video, u_controls_8_2)
    ee_positions_4_controlled_2, qpos_4_controlled_2, _ = simulate(model_4, data_4_inp, duration_controlled_2, ee_site_id_4, final_qpos_4, final_qvel_4, create_video, u_controls_4_2)
    ee_positions_2_controlled_2, qpos_2_controlled_2, _ = simulate(model_2, data_2_inp, duration_controlled_2, ee_site_id_2, final_qpos_2, final_qvel_2, create_video, u_controls_2_2)

    ee_positions_16_controlled.extend(ee_positions_16_controlled_2)
    ee_positions_8_controlled.extend(ee_positions_8_controlled_2)
    ee_positions_4_controlled.extend(ee_positions_4_controlled_2)
    ee_positions_2_controlled.extend(ee_positions_2_controlled_2)

    qpos_16_controlled.extend(qpos_16_controlled_2)
    qpos_8_controlled.extend(qpos_8_controlled_2)
    qpos_4_controlled.extend(qpos_4_controlled_2)
    qpos_2_controlled.extend(qpos_2_controlled_2)

    final_qpos_16 = qpos_16_controlled[-1]
    final_qpos_8 = qpos_8_controlled[-1]
    final_qpos_4 = qpos_4_controlled[-1]
    final_qpos_2 = qpos_2_controlled[-1]

    duration_controlled += duration_controlled_2

    plot_chains([link_length_16, link_length_8, link_length_4, link_length_2], [final_qpos_16, final_qpos_8, final_qpos_4, final_qpos_2], fontsizes)
    
# Optionally, include a title or save the figure
plt.title("Final Configuration of Each Chain")
#plt.show()

# Extract x and z positions over time for each controlled system
time_controlled = np.linspace(0, duration_controlled, len(ee_positions_16_controlled))

ee_positions_16_controlled_x = [pos[0] for pos in ee_positions_16_controlled]
ee_positions_16_controlled_z = [pos[2] for pos in ee_positions_16_controlled]

ee_positions_8_controlled_x = [pos[0] for pos in ee_positions_8_controlled]
ee_positions_8_controlled_z = [pos[2] for pos in ee_positions_8_controlled]

ee_positions_4_controlled_x = [pos[0] for pos in ee_positions_4_controlled]
ee_positions_4_controlled_z = [pos[2] for pos in ee_positions_4_controlled]

ee_positions_2_controlled_x = [pos[0] for pos in ee_positions_2_controlled]
ee_positions_2_controlled_z = [pos[2] for pos in ee_positions_2_controlled]

# Plot endeffector x and z positions over time in a 2x1 subplot
fig, axes = plt.subplots(2, 1, figsize=(10, 8))
axes[0].plot(time_controlled, ee_positions_16_controlled_x, label="16-link (controlled)")
axes[0].plot(time_controlled, ee_positions_8_controlled_x, label="8-link (controlled)")
axes[0].plot(time_controlled, ee_positions_4_controlled_x, label="4-link (controlled)")
axes[0].plot(time_controlled, ee_positions_2_controlled_x, label="2-link (controlled)")
axes[0].set_title('Controlled Endeffector X Position Over Time')
axes[0].set_xlabel('Time (s)')
axes[0].set_ylabel('X Position')
axes[0].legend()

axes[1].plot(time_controlled, ee_positions_16_controlled_z, label="16-link (controlled)")
axes[1].plot(time_controlled, ee_positions_8_controlled_z, label="8-link (controlled)")
axes[1].plot(time_controlled, ee_positions_4_controlled_z, label="4-link (controlled)")
axes[1].plot(time_controlled, ee_positions_2_controlled_z, label="2-link (controlled)")
axes[1].set_title('Controlled Endeffector Z Position Over Time')
axes[1].set_xlabel('Time (s)')
axes[1].set_ylabel('Z Position')
axes[1].legend()

plt.tight_layout()
plt.show()

# Plot development in the x-z plane for controlled systems
plt.figure(figsize=(10, 6))
plt.plot(ee_positions_16_controlled_x, ee_positions_16_controlled_z, label="16-link (controlled)")
plt.plot(ee_positions_8_controlled_x, ee_positions_8_controlled_z, label="8-link (controlled)")
plt.plot(ee_positions_4_controlled_x, ee_positions_4_controlled_z, label="4-link (controlled)")
plt.plot(ee_positions_2_controlled_x, ee_positions_2_controlled_z, label="2-link (controlled)")
plt.title('Controlled Endeffector Position in X-Z Plane')
plt.xlabel('X Position')
plt.ylabel('Z Position')
plt.legend()
plt.grid()
plt.show()

# Check behaviour under sinusoidal control inputs
sinusoidal_control = True
if sinusoidal_control:
    # Define sinusoidal control inputs for u_16_1 and u_16_2
    frequency = 1.0  # Frequency of sinusoidal inputs (Hz)
    amplitude = 0.4  # Amplitude of sinusoidal inputs
    duration_sinusoidal = 6.0  # Duration of the simulation (seconds)
    time_sinusoidal = np.linspace(0, duration_sinusoidal, int(np.ceil(duration_sinusoidal/model_16.opt.timestep)))  # Time array for sinusoidal inputs

    u_16_1_sinusoidal = amplitude * np.sin(2 * np.pi * frequency * time_sinusoidal)
    u_16_2_sinusoidal = amplitude * np.cos(2 * np.pi * frequency * time_sinusoidal)

    # compute controls
    u_controls_2 = []
    u_controls_4 = []
    u_controls_8 = []
    u_controls_16 = []
    for u_1, u_2 in zip(u_16_1_sinusoidal, u_16_2_sinusoidal):
        u_controls_2.append(compute_controls(u_1, u_2, link_length_16, link_length_2, stiffness_16, stiffness_2, gear, 16, 2))
        u_controls_4.append(compute_controls(u_1, u_2, link_length_16, link_length_4, stiffness_16, stiffness_4, gear, 16, 4))
        u_controls_8.append(compute_controls(u_1, u_2, link_length_16, link_length_8, stiffness_16, stiffness_8, gear, 16, 8))
        u_controls_16.append(np.array([u_1, u_2]))
    
    if not map_controls:
        u_controls_2 = u_controls_16
        u_controls_4 = u_controls_16
        u_controls_8 = u_controls_16
    # simulate    
    ee_positions_16_sinusoidal, q_positions_16_sinusoidal, _ = simulate(model_16, deepcopy(data_16), duration_sinusoidal, ee_site_id_16, qpos_16_ss, qvel_16_init, False, u_controls_16)
    ee_positions_8_sinusoidal, q_positions_8_sinusoidal, _ = simulate(model_8, deepcopy(data_8), duration_sinusoidal, ee_site_id_8, qpos_8_ss, qvel_8_init, False, u_controls_8)
    ee_positions_4_sinusoidal, q_positions_4_sinusoidal, _ = simulate(model_4, deepcopy(data_4), duration_sinusoidal, ee_site_id_4, qpos_4_ss, qvel_4_init, False, u_controls_4)
    ee_positions_2_sinusoidal, q_positions_2_sinusoidal, _ = simulate(model_2, deepcopy(data_2), duration_sinusoidal, ee_site_id_2, qpos_2_ss, qvel_2_init, False, u_controls_2)

    # Extract x and z positions
    ee_16_x_sinusoidal = [pos[0] for pos in ee_positions_16_sinusoidal]
    ee_16_z_sinusoidal = [pos[2] for pos in ee_positions_16_sinusoidal]

    ee_8_x_sinusoidal = [pos[0] for pos in ee_positions_8_sinusoidal]
    ee_8_z_sinusoidal = [pos[2] for pos in ee_positions_8_sinusoidal]

    ee_4_x_sinusoidal = [pos[0] for pos in ee_positions_4_sinusoidal]
    ee_4_z_sinusoidal = [pos[2] for pos in ee_positions_4_sinusoidal]

    ee_2_x_sinusoidal = [pos[0] for pos in ee_positions_2_sinusoidal]
    ee_2_z_sinusoidal = [pos[2] for pos in ee_positions_2_sinusoidal]

    # Plot the sinusoidal trajectories in the X-Z plane
    plt.figure(figsize=(10, 6))
    plt.plot(ee_16_x_sinusoidal, ee_16_z_sinusoidal, label="16-link (sinusoidal)")
    plt.plot(ee_8_x_sinusoidal, ee_8_z_sinusoidal, label="8-link (sinusoidal)")
    plt.plot(ee_4_x_sinusoidal, ee_4_z_sinusoidal, label="4-link (sinusoidal)")
    plt.plot(ee_2_x_sinusoidal, ee_2_z_sinusoidal, label="2-link (sinusoidal)")
    plt.title('Sinusoidal Endeffector Position in X-Z Plane')
    plt.xlabel('X Position')
    plt.ylabel('Z Position')
    plt.legend()
    plt.grid()
    plt.show()

    # Plot X and Z positions over time
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    # X positions over time
    axes[0].plot(time_sinusoidal, ee_16_x_sinusoidal, label="16-link (sinusoidal)")
    axes[0].plot(time_sinusoidal, ee_8_x_sinusoidal, label="8-link (sinusoidal)")
    axes[0].plot(time_sinusoidal, ee_4_x_sinusoidal, label="4-link (sinusoidal)")
    axes[0].plot(time_sinusoidal, ee_2_x_sinusoidal, label="2-link (sinusoidal)")
    axes[0].set_title('Endeffector X Position Over Time (Sinusoidal Controls)')
    axes[0].set_xlabel('Time (s)')
    axes[0].set_ylabel('X Position')
    axes[0].legend()
    axes[0].grid()

    # Z positions over time
    axes[1].plot(time_sinusoidal, ee_16_z_sinusoidal, label="16-link (sinusoidal)")
    axes[1].plot(time_sinusoidal, ee_8_z_sinusoidal, label="8-link (sinusoidal)")
    axes[1].plot(time_sinusoidal, ee_4_z_sinusoidal, label="4-link (sinusoidal)")
    axes[1].plot(time_sinusoidal, ee_2_z_sinusoidal, label="2-link (sinusoidal)")
    axes[1].set_title('Endeffector Z Position Over Time (Sinusoidal Controls)')
    axes[1].set_xlabel('Time (s)')
    axes[1].set_ylabel('Z Position')
    axes[1].legend()
    axes[1].grid()

    plt.tight_layout()
    plt.show()

    # create plot with both autonomous and sinusoidal controlled trajectories for visualization
    end_index_autonomous = 500
    end_index_controlled = 1000

    # Define the number of subplots per row
    ncols = 2

    # Create the main figure
    fig = plt.figure(figsize=(12, 10), constrained_layout=True)

    # Create two subfigures, one for each row
    subfigs = fig.subfigures(nrows=2, ncols=1)

    # First row: Autonomous Trajectories
    subfigs[0].suptitle('Autonomous Trajectories', fontsize=22)
    axs_top = subfigs[0].subplots(1, ncols)

    # Plot for X position (Autonomous)
    axs_top[0].plot(time[:end_index_autonomous], ee_positions_16_x[:end_index_autonomous], label="16 links")
    axs_top[0].plot(time[:end_index_autonomous], ee_positions_8_x[:end_index_autonomous], label="8 links")
    axs_top[0].plot(time[:end_index_autonomous], ee_positions_4_x[:end_index_autonomous], label="4 links")
    axs_top[0].plot(time[:end_index_autonomous], ee_positions_2_x[:end_index_autonomous], label="2 links")
    axs_top[0].set_xlabel('Time (s)')
    axs_top[0].set_ylabel('x')
    axs_top[0].grid(True)
    axs_top[0].legend()

    # Plot for Z position (Autonomous)
    axs_top[1].plot(time[:end_index_autonomous], ee_positions_16_z[:end_index_autonomous])
    axs_top[1].plot(time[:end_index_autonomous], ee_positions_8_z[:end_index_autonomous])
    axs_top[1].plot(time[:end_index_autonomous], ee_positions_4_z[:end_index_autonomous])
    axs_top[1].plot(time[:end_index_autonomous], ee_positions_2_z[:end_index_autonomous])
    axs_top[1].set_xlabel('Time (s)')
    axs_top[1].set_ylabel('z')
    axs_top[1].grid(True)

    # Second row: Controlled Trajectories under Sinusoidal u
    subfigs[1].suptitle('Controlled Trajectories under Sinusoidal $\\mathbf{u}$', fontsize=22)
    axs_bottom = subfigs[1].subplots(1, ncols)

    # Plot for X position (Controlled)
    axs_bottom[0].plot(time_sinusoidal[:end_index_controlled], ee_16_x_sinusoidal[:end_index_controlled])
    axs_bottom[0].plot(time_sinusoidal[:end_index_controlled], ee_8_x_sinusoidal[:end_index_controlled])
    axs_bottom[0].plot(time_sinusoidal[:end_index_controlled], ee_4_x_sinusoidal[:end_index_controlled])
    axs_bottom[0].plot(time_sinusoidal[:end_index_controlled], ee_2_x_sinusoidal[:end_index_controlled])
    axs_bottom[0].set_xlabel('Time (s)')
    axs_bottom[0].set_ylabel('x')
    axs_bottom[0].grid(True)

    # Plot for Z position (Controlled)
    axs_bottom[1].plot(time_sinusoidal[:end_index_controlled], ee_16_z_sinusoidal[:end_index_controlled])
    axs_bottom[1].plot(time_sinusoidal[:end_index_controlled], ee_8_z_sinusoidal[:end_index_controlled])
    axs_bottom[1].plot(time_sinusoidal[:end_index_controlled], ee_4_z_sinusoidal[:end_index_controlled])
    axs_bottom[1].plot(time_sinusoidal[:end_index_controlled], ee_2_z_sinusoidal[:end_index_controlled])
    axs_bottom[1].set_xlabel('Time (s)')
    axs_bottom[1].set_ylabel('z')
    axs_bottom[1].grid(True)

    # Display the plot
    os.makedirs("plots", exist_ok=True)
    filename = f"plots/autonomous_and_controlled_trajectories.pgf"
    plt.savefig(filename, bbox_inches="tight")
    plt.show()

    ani = animate_chains(
        q_trajectories=[q_positions_16_sinusoidal, q_positions_8_sinusoidal, q_positions_4_sinusoidal, q_positions_2_sinusoidal],
        link_lengths=[link_length_16, link_length_8, link_length_4, link_length_2],
        fontsizes=fontsizes,
        duration=duration
    )

# --- simulate mixed trajectories with sinusoidal controls, where we map an intermediate state ---
# Symbolic setup for velocity mapping using CasADi and chain rule
q_high_dim_16 = ca.SX.sym('q_high_dim_16', 16)  # Symbolic q_high for 16-link
q_high_dot_16 = ca.SX.sym('q_high_dot_16', 16)  # Symbolic q_high_dot for 16-link

q_8 = compute_q_casadi(q_high_dim_16, link_length_16, link_length_8, 8)
jacobian_q_8_q_16 = ca.jacobian(q_8, q_high_dim_16)  # Jacobian of q_8 w.r.t. q_16
q_dot_8 = jacobian_q_8_q_16 @ q_high_dot_16  # Velocity mapping for 8-link
q_low_dot_fn_8 = ca.Function('q_low_dot_fn_8', [q_high_dim_16, q_high_dot_16], [q_dot_8])  # CasADi function

q_high_dot_8 = ca.SX.sym('q_high_dot_8', 8)  # Symbolic q_high_dot for 8-link

q_4 = compute_q_casadi(q_high_dim_16, link_length_16, link_length_4, 4)
jacobian_q_4_q_16 = ca.jacobian(q_4, q_high_dim_16)  # Jacobian of q_4 w.r.t. q_16
q_dot_4 = jacobian_q_4_q_16 @ q_high_dot_16  # Velocity mapping for 4-link
q_low_dot_fn_4 = ca.Function('q_low_dot_fn_4', [q_high_dim_16, q_high_dot_16], [q_dot_4])  # CasADi function

q_high_dot_4 = ca.SX.sym('q_high_dot_4', 4)  # Symbolic q_high_dot for 4-link

q_2 = compute_q_casadi(q_high_dim_16, link_length_4, link_length_2, 2)
jacobian_q_2_q_16 = ca.jacobian(q_2, q_high_dim_16)  # Jacobian of q_2 w.r.t. q_4
q_dot_2 = jacobian_q_2_q_16 @ q_high_dot_16  # Velocity mapping for 2-link
q_low_dot_fn_2 = ca.Function('q_low_dot_fn_2', [q_high_dim_16, q_high_dot_16], [q_dot_2])  # CasADi function

# Define sinusoidal control inputs for the 16-link chain
frequency = 1.0  # Frequency of the sinusoidal input (Hz)
amplitude = 0.47 # Amplitude of the sinusoidal input
time_array = np.linspace(0, duration_controlled, int(np.ceil(duration_controlled/model_16.opt.timestep)))  # Time array for the controls

u_16_1_sinusoidal = amplitude * np.sin(2 * np.pi * frequency * time_array)
u_16_2_sinusoidal = amplitude * np.cos(2 * np.pi * frequency * time_array)

u_controls_16 = np.vstack((u_16_1_sinusoidal, u_16_2_sinusoidal)).T

ee_positions_16, qpos_16_controlled, qvel_16_controlled = simulate(
    model_16, data_16, duration_controlled, ee_site_id_16, qpos_16_ss, qvel_16_init, create_video=False, u_controls=u_controls_16
)
    
# Extract the state in the middle of the sinusoidal trajectory
middle_index = len(qpos_16_controlled) // 2 -1
qpos_16_mid = qpos_16_controlled[middle_index]
qvel_16_mid = qvel_16_controlled[middle_index]

# Map the middle state to lower-dimensional systems using compute_q and CasADi
qpos_8_mid = compute_q(qpos_16_mid, link_length_16, link_length_8, 8)
qpos_4_mid = compute_q(qpos_16_mid, link_length_16, link_length_4, 4)
qpos_2_mid = compute_q(qpos_16_mid, link_length_16, link_length_2, 2)

# Compute velocities for the 8, 4, and 2-link chains
qvel_8_mid = q_low_dot_fn_8(qpos_16_mid, qvel_16_mid)
qvel_4_mid = q_low_dot_fn_4(qpos_16_mid, qvel_16_mid)
qvel_2_mid = q_low_dot_fn_2(qpos_16_mid, qvel_16_mid)

# compute_controls 
if map_controls:
    u_controls_8 = []
    u_controls_4 = []
    u_controls_2 = []
    for u_1, u_2 in zip(u_16_1_sinusoidal[middle_index+1:], u_16_2_sinusoidal[middle_index+1:]):
        # Map sinusoidal controls to lower-dimensional systems
        u_controls_8.append(compute_controls(u_1, u_2, link_length_16, link_length_8, stiffness_16, stiffness_8, gear, 16, 8))
        u_controls_4.append(compute_controls(u_1, u_2, link_length_16, link_length_4, stiffness_16, stiffness_4, gear, 16, 4))
        u_controls_2.append(compute_controls(u_1, u_2, link_length_16, link_length_2, stiffness_16, stiffness_2, gear, 16, 2))
else:
    u_controls_8 = u_controls_16[middle_index+1:]
    u_controls_4 = u_controls_16[middle_index+1:]
    u_controls_2 = u_controls_16[middle_index+1:]

# simulate
ee_positions_8, q_pos_8_switch, _ = simulate(
    model_8, data_8, duration_controlled/2, ee_site_id_8, qpos_8_mid, np.array(qvel_8_mid).squeeze(), create_video=False, u_controls=u_controls_8
)
ee_positions_4, q_pos_4_switch, _ = simulate(
    model_4, data_4, duration_controlled/2, ee_site_id_4, qpos_4_mid, np.array(qvel_4_mid).squeeze(), create_video=False, u_controls=u_controls_4
)
ee_positions_2, q_pos_2_switch, _ = simulate(
    model_2, data_2, duration_controlled/2, ee_site_id_2, qpos_2_mid, np.array(qvel_2_mid).squeeze(), create_video=False, u_controls=u_controls_2
)

# Plot the resulting trajectories in the X-Z plane
plt.figure(figsize=(10, 6))
plt.plot([pos[0] for pos in ee_positions_8], [pos[2] for pos in ee_positions_8], label="8-link (mapped)")
plt.plot([pos[0] for pos in ee_positions_4], [pos[2] for pos in ee_positions_4], label="4-link (mapped)")
plt.plot([pos[0] for pos in ee_positions_2], [pos[2] for pos in ee_positions_2], label="2-link (mapped)")
plt.title('Mapped Endeffector Position in X-Z Plane with Sinusoidal Control')
plt.xlabel('X Position')
plt.ylabel('Z Position')
plt.legend()
plt.grid()
plt.show()

# Plot the X and Z positions over time
time_total = np.linspace(-duration_controlled/2, duration_controlled/2, len(ee_positions_16))
time_mapped = np.linspace(0, duration_controlled/2, len(ee_positions_8))

fig, axes = plt.subplots(2, 1, figsize=(10, 8))
# X positions
axes[0].plot(time_total, [pos[0] for pos in ee_positions_16], label="16-link (mapped)")
axes[0].plot(time_mapped, [pos[0] for pos in ee_positions_8], label="8-link (mapped)")
axes[0].plot(time_mapped, [pos[0] for pos in ee_positions_4], label="4-link (mapped)")
axes[0].plot(time_mapped, [pos[0] for pos in ee_positions_2], label="2-link (mapped)")
axes[0].set_title('Endeffector X Position Over Time (Mapped, Sinusoidal)')
axes[0].set_xlabel('Time (s)')
axes[0].set_ylabel('X Position')
axes[0].legend()
axes[0].grid()

# Z positions
axes[1].plot(time_total, [pos[2] for pos in ee_positions_16], label="16-link (mapped)")
axes[1].plot(time_mapped, [pos[2] for pos in ee_positions_8], label="8-link (mapped)")
axes[1].plot(time_mapped, [pos[2] for pos in ee_positions_4], label="4-link (mapped)")
axes[1].plot(time_mapped, [pos[2] for pos in ee_positions_2], label="2-link (mapped)")
axes[1].set_title('Endeffector Z Position Over Time (Mapped, Sinusoidal)')
axes[1].set_xlabel('Time (s)')
axes[1].set_ylabel('Z Position')
axes[1].legend()
axes[1].grid()

plt.tight_layout()
plt.show()

ani = animate_chains(
    q_trajectories=[qpos_16_controlled[middle_index+1:], q_pos_8_switch, q_pos_4_switch, q_pos_2_switch],
    link_lengths=[link_length_16, link_length_8, link_length_4, link_length_2],
    fontsizes=fontsizes,
    duration=duration
)
