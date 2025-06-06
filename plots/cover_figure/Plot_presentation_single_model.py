"""
This file creates a figure used in the motivation of the thesis presentation.
"""
import matplotlib.pyplot as plt
import numpy as np
from acados_template import latexify_plot

# Latex-style plotting
latexify_plot(fontsize=15)

# === Signal setup ===
t_total = 7.0
num_points = 1000

t = np.linspace(0, t_total, num_points)

# Slowly varying component (same as original)
slow_oscillation = (
    0.6 * np.sin(0.4 * np.pi * t + 0.3) +
    0.15 * np.sin(0.9 * np.pi * t)
)

# Fast-varying component (same as original)
fast_oscillation = (
    0.2 * np.sin(7 * np.pi * t + 0.1) +
    0.1 * np.sin(11 * np.pi * t + np.pi / 5) +
    0.05 * np.sin(17 * np.pi * t + 1.2)
)

# Noise component
tmp_rng = np.random.default_rng(seed=42)
noise = 0.02 * np.convolve(tmp_rng.normal(0, 1, len(t)), np.ones(10)/10, mode='same')

# Combine into full, fast-varying signal (old 'first phase')
offset = 1.3
signal = offset + slow_oscillation + fast_oscillation + noise

# === Constant-step Riemann integration ===
n_steps = 90
dt = t_total / n_steps
t_riemann = np.linspace(0, t_total, n_steps + 1)

t_left = t_riemann[:-1]
signal_left = np.interp(t_left, t, signal)
areas = signal_left * dt
riemann_sum = np.sum(areas)
print(f"Approximate integral (Left Riemann sum) = {riemann_sum:.4f}")

# === Plotting ===
fig, ax = plt.subplots()

# Axes arrows
head_len = 0.3
shaft_len = 2.7  # arrow shaft length for y-axis
ax.arrow(-0.3, 0, t_total + 1.0, 0,
         head_width=head_len/2, head_length=head_len, fc='black', ec='black')
ax.arrow(-0.3, 0, 0, shaft_len,
         head_width=head_len/2, head_length=head_len, fc='black', ec='black')

# Axis labels
ax.text(t_total + 1.0, -0.2, r'$t$', va='top', ha='center', fontsize=20)
ax.text(-0.6, shaft_len, r'$\dot{x}$', va='bottom', ha='center', fontsize=20)

# Tick marks at t0 and tN
tick_height = 0.1
ax.plot([0, 0], [-tick_height, tick_height], color='black', linewidth=1.0)
ax.text(0, -0.2, r'$t_0$', va='top', ha='center', fontsize=20)
ax.plot([t_total, t_total], [-tick_height, tick_height], color='black', linewidth=1.0)
ax.text(t_total, -0.2, r'$t_{N}$', va='top', ha='center', fontsize=20)

# Dashed vertical lines
dashed_line_height = 3.0
ax.plot([0, 0], [0, dashed_line_height], linestyle='--', color='black', linewidth=1.0)
ax.plot([t_total, t_total], [0, dashed_line_height], linestyle='--', color='black', linewidth=1.0)

# Signal curve
ax.plot(t, signal, lw=1.8, color='black')

# Riemann rectangles
for t0, s0 in zip(t_left, signal_left):
    ax.bar(t0, s0, width=dt, align='edge', color='orange', alpha=0.5,
           edgecolor='black', linewidth=0.3)

# Final formatting
ax.set_xlim(-0.6, t_total + 1.5)
ax.set_ylim(-0.6, dashed_line_height + 0.5)
ax.set_aspect('equal')
ax.axis('off')

plt.show()
