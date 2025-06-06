"""
This file creates the cover figure for the fading-fidelity MPC approach for slow-fast systems.
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.optimize import root_scalar
from acados_template import latexify_plot

latexify_plot(fontsize=15)

# === Helper: compute geometric schedule ===
def compute_geometric_schedule(dt0, n_steps, total_length):
    def geometric_sum_error(r):
        if r == 1.0:
            return dt0 * n_steps - total_length
        return dt0 * (1 - r**n_steps) / (1 - r) - total_length

    result = root_scalar(geometric_sum_error, bracket=[1.00001, 100], method='brentq')
    if not result.converged:
        raise RuntimeError("Could not find suitable growth rate.")
    r = result.root
    steps = dt0 * r**np.arange(n_steps)
    t_riemann = np.cumsum(np.concatenate(([0.0], steps)))
    return t_riemann

# === Signal setup ===
t_arrow_length = 7.0
x_dot_arrow_length = 3.0
dashed_line_height = 3.0  # <-- adjustable height for all dashed lines

t = np.linspace(0, t_arrow_length, 1000)

slow_oscillation = 0.6 * np.sin(0.4 * np.pi * t + 0.3) + 0.15 * np.sin(0.9 * np.pi * t)
fast_oscillation = (
    0.2 * np.sin(7 * np.pi * t + 0.1) +
    0.1 * np.sin(11 * np.pi * t + np.pi / 5) +
    0.05 * np.sin(17 * np.pi * t + 1.2)
)
rng = np.random.default_rng(seed=42)
noise = 0.02 * np.convolve(rng.normal(0, 1, len(t)), np.ones(10)/10, mode='same')

# === Riemann integration ===
start, end = t[0], t[-1]
dt0 = 0.05
n_steps = 30

t_riemann = compute_geometric_schedule(dt0, n_steps, end - start)
t_riemann += start
t_riemann = t_riemann[t_riemann <= end]

switch_step = 18
t_switch = t_riemann[switch_step]

# Smooth transition
transition_width = 0.03
transition = 0.5 * (1 - np.tanh((t - t_switch) / transition_width))

# Signal
signal = 1.3 + slow_oscillation + (fast_oscillation + noise) * transition

# Riemann sum (left)
t_leftpoints = t_riemann[:-1]
signal_leftpoints = np.interp(t_leftpoints, t, signal)
areas = signal_leftpoints * np.diff(t_riemann)
riemann_sum = np.sum(areas)

print(f"Approximate integral (Left Riemann sum) = {riemann_sum:.4f}")

# === Plotting ===
fig, ax = plt.subplots()

# Axes
head_len = 0.3
shaft_len = x_dot_arrow_length - head_len
ax.arrow(-0.3, 0, t_arrow_length + 1.0, 0, head_width=head_len/2, head_length=head_len, fc='black', ec='black')
ax.arrow(-0.3, 0, 0, shaft_len, head_width=head_len/2, head_length=head_len, fc='black', ec='black')

# Axis labels
ax.text(t_arrow_length + 1.0, -0.2, r'$t$', va='top', ha='center', fontsize=20)
ax.text(-0.6, x_dot_arrow_length, r'$\dot{x}$', va='bottom', ha='center', fontsize=20)

# Markers at t_0, t_switch, and T
tick_height = 0.1
ax.plot([0, 0], [-tick_height, tick_height], color='black', linewidth=1.0)
ax.text(0, -0.2, r'$t_0$', va='top', ha='center', fontsize=20)

ax.plot([t_switch, t_switch], [-tick_height, tick_height], color='black', linewidth=1.0)
ax.text(t_switch, -0.2, r'$t_{\mathrm{switch}}$', va='top', ha='center', fontsize=20)

ax.plot([t[-1], t[-1]], [-tick_height, tick_height], color='black', linewidth=1.0)
ax.text(t[-1], -0.2, r'$t_{N}$', va='top', ha='center', fontsize=20)

# Dashed vertical lines at t_0, t_switch, and t_N
ax.plot([0, 0], [0, dashed_line_height], linestyle='--', color='black', linewidth=1.0)
ax.plot([t_switch, t_switch], [0, dashed_line_height], linestyle='--', color='black', linewidth=1.0)
ax.plot([t[-1], t[-1]], [0, dashed_line_height], linestyle='--', color='black', linewidth=1.0)

# Signal curve
ax.plot(t, signal, lw=1.8, color='black')

# Riemann rectangles
for i in range(len(t_leftpoints)):
    color = 'orange' if i < switch_step else 'lightgreen'
    ax.bar(t_leftpoints[i], signal_leftpoints[i], width=np.diff(t_riemann)[i],
           align='edge', color=color, alpha=0.5, edgecolor='black', linewidth=0.3)

# Final formatting
ax.set_xlim(-0.6, t_arrow_length + 1.5)
ax.set_ylim(-0.6, dashed_line_height + 0.5)
ax.set_aspect('equal')
ax.axis('off')

plt.show()

plt.show()


