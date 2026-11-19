# Mohammad Zare
# Slider-crank dynamic balancing (optional bonus project)
#
# Computes the shaking force of a slider-crank mechanism (rotating part plus
# primary and secondary reciprocating parts), then adds a crank counterweight
# that balances all of the rotating mass and a fraction x of the reciprocating
# mass. It searches x from 0 to 1 and picks the value that minimizes the
# shaking force over one crank revolution.

import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, solve

# Define parameters (SI units: kg, m)
connecting_rod_mass = 15.4
crank_mass = 4.54
piston_mass = 9.1
connecting_rod = 0.356          # length l of the connecting rod
stroke = 0.204
crank = stroke / 2              # crank radius r
cm_crank = 0.0381               # distance of the crank's center of mass from the main bearing
cm_connecting_rod = 0.102       # distance of the rod's center of mass from the crank pin
w2_rpm = 1000
w2 = w2_rpm * (2 * np.pi / 60)  # crank speed in rad/s
theta_values = np.linspace(0, 2 * np.pi, 1000)


def crank_mass_func(cm_crank, crank_mass, crank):
    """Equivalent lumped mass of the crank, referred to the crank pin."""
    return (cm_crank * crank_mass) / crank


def connecting_rod_masses(connecting_rod_mass, cm_connecting_rod, connecting_rod):
    """Split the connecting rod's mass into two lumped masses at its two
    ends (crank pin side and piston pin side) that preserve the rod's
    total mass and center of mass location."""
    m2, m3 = symbols('m2 m3')
    eq1 = m2 + m3 - connecting_rod_mass
    eq2 = (m2 * cm_connecting_rod) - (m3 * (connecting_rod - cm_connecting_rod))
    solution = solve((eq1, eq2), (m2, m3))
    # Convert SymPy numbers to plain floats so NumPy does fast numeric math later
    return float(solution[m2]), float(solution[m3])


# Calculate lumped masses
crank_mass1 = crank_mass_func(cm_crank, crank_mass, crank)
connecting_rod_mass1, connecting_rod_mass2 = connecting_rod_masses(
    connecting_rod_mass, cm_connecting_rod, connecting_rod
)
rotating_mass = crank_mass1 + connecting_rod_mass1        # m_A, moves on a circle at the crank pin
reciprocating_mass = piston_mass + connecting_rod_mass2   # m_B, moves back and forth with the piston


def forces_function(theta_values):
    """Unbalanced shaking force (no counterweight) for every crank angle.
    Complex number: real part = horizontal, imaginary part = vertical."""
    rotating = rotating_mass * crank * w2 ** 2 * np.exp(1j * theta_values)
    primary = reciprocating_mass * crank * w2 ** 2 * np.cos(theta_values)
    secondary = reciprocating_mass * crank * w2 ** 2 * (crank / connecting_rod) * np.cos(2 * theta_values)
    return rotating + primary + secondary


def shaking_force_function(x, theta_values):
    """Shaking force with a counterweight opposite the crank pin that balances
    all of the rotating mass plus a fraction x (0 to 1) of the reciprocating mass."""
    counterweight = -(rotating_mass + x * reciprocating_mass) * crank * w2 ** 2 * np.exp(1j * theta_values)
    return forces_function(theta_values) + counterweight


def mean_force_magnitude(force_values, theta_values):
    """Average shaking-force magnitude over one revolution
    (trapezoidal integral of |F| over theta, divided by 2*pi)."""
    force_magnitudes = np.abs(force_values)
    delta_theta = np.diff(theta_values)
    integral = np.sum(0.5 * (force_magnitudes[:-1] + force_magnitudes[1:]) * delta_theta)
    return integral / (2 * np.pi)


# Baseline (unbalanced) shaking force
force_values = forces_function(theta_values)
original_mean = mean_force_magnitude(force_values, theta_values)

# Plot the baseline force locus
plt.figure(figsize=(8, 8))
plt.plot(np.real(force_values), np.imag(force_values), label='Unbalanced shaking force')
plt.title('Horizontal and vertical parts of the shaking force (no counterweight)')
plt.xlabel('Horizontal force (N)')
plt.ylabel('Vertical force (N)')
plt.axis('equal')
plt.legend()
plt.grid(True)
plt.show()

# Search balancing factors x from 0 to 1 and keep the one with the smallest mean shaking force
x_values = np.linspace(0, 1, 101)
mean_forces = [mean_force_magnitude(shaking_force_function(x, theta_values), theta_values) for x in x_values]
best_index = int(np.argmin(mean_forces))
final_x = x_values[best_index]
final_fs = shaking_force_function(final_x, theta_values)

# Plot how the mean shaking force depends on x
plt.figure(figsize=(8, 5))
plt.plot(x_values, mean_forces)
plt.axvline(final_x, color='red', linestyle='dashed', label=f'best x = {final_x:.2f}')
plt.title('Mean shaking force vs. balancing factor x')
plt.xlabel('Balancing factor x (fraction of reciprocating mass balanced)')
plt.ylabel('Mean shaking force (N)')
plt.legend()
plt.grid(True)
plt.show()

# Plot original vs. balanced shaking force
plt.figure(figsize=(10, 8))
plt.plot(np.real(force_values), np.imag(force_values), label='Unbalanced', color='blue')
plt.plot(np.real(final_fs), np.imag(final_fs), label=f'Balanced, x = {final_x:.2f}', linestyle='dashed', color='red')
plt.title('Shaking force before and after balancing')
plt.xlabel('Horizontal force (N)')
plt.ylabel('Vertical force (N)')
plt.axis('equal')
plt.legend()
plt.grid(True)
plt.show()

print(f"Best balancing factor x = {final_x:.2f}")
print(f"Mean shaking force: {original_mean:.0f} N unbalanced -> {mean_forces[best_index]:.0f} N balanced "
      f"({100 * (1 - mean_forces[best_index] / original_mean):.0f}% lower)")
print(f"Peak shaking force: {np.abs(force_values).max():.0f} N unbalanced -> {np.abs(final_fs).max():.0f} N balanced")
