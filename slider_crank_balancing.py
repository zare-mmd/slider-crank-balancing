# Mohammad Zare
# Slider-crank dynamic balancing (optional bonus project)
#
# Computes the primary and secondary shaking forces of a slider-crank
# mechanism, then searches over a range of counterweight balancing
# factors (x) to find the value that minimizes the shaking force,
# measured as the area enclosed by the force-ellipse curve.

import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, solve

# Define parameters
connecting_rod_mass = 15.4
crank_mass = 4.54
piston_mass = 9.1
connecting_rod = 0.356
stroke = 0.204
crank = stroke / 2
cm_crank = 0.0381
cm_connecting_rod = 0.102
w2_rpm = 1000
w2 = w2_rpm * (2 * np.pi / 60)
theta_values = np.linspace(0, 2 * np.pi, 1000)


def crank_mass_func(cm_crank, crank_mass, crank):
    """Equivalent lumped mass of the crank, referred to the crank pin."""
    return (cm_crank * crank_mass) / crank


def connecting_rod_masses(connecting_rod_mass, cm_connecting_rod, connecting_rod):
    """Split the connecting rod's mass into two lumped masses at its two
    ends (crank pin side and piston pin side) that preserve the rod's
    total mass and center of mass location."""
    m2, m3, theta = symbols('m2 m3 theta')
    eq1 = m2 + m3 - connecting_rod_mass
    eq2 = (m2 * cm_connecting_rod) - (m3 * (connecting_rod - cm_connecting_rod))
    solution = solve((eq1, eq2), (m2, m3))
    return solution[m2], solution[m3]


# Calculate lumped masses
crank_mass1 = crank_mass_func(cm_crank, crank_mass, crank)
connecting_rod_mass1, connecting_rod_mass2 = connecting_rod_masses(
    connecting_rod_mass, cm_connecting_rod, connecting_rod
)


def forces_function(piston_mass, connecting_rod_mass2, crank, w2, connecting_rod,
                     crank_mass1, connecting_rod_mass1, theta_values):
    """Total shaking force (rotating + primary + secondary components),
    represented as a complex number per crank angle so the horizontal
    and vertical components can be plotted directly."""
    primary = ((piston_mass + connecting_rod_mass2) * crank * (w2 ** 2)) * np.cos(theta_values)
    secondary = ((piston_mass + connecting_rod_mass2) * (crank ** 2) * (w2 ** 2) / connecting_rod) * np.cos(2 * theta_values)
    fc = (crank_mass1 + connecting_rod_mass1) * crank * (w2 ** 2)
    forces = fc * np.exp(1j * theta_values) + primary + secondary
    return np.array([complex(force.evalf()) for force in forces], dtype=np.complex128)


# Calculate baseline (unbalanced) force values
force_values = forces_function(
    piston_mass, connecting_rod_mass2, crank, w2, connecting_rod,
    crank_mass1, connecting_rod_mass1, theta_values
)

# Plot the baseline force locus
plt.figure(figsize=(8, 8))
plt.plot(np.real(force_values), np.imag(force_values), label='Forces')
plt.title('Horizontal and Vertical Parts of Forces')
plt.xlabel('Horizontal Forces')
plt.ylabel('Vertical Forces')
plt.legend()
plt.grid(True)
plt.show()


def shaking_force_function(x, piston_mass, connecting_rod_mass2, crank, w2, connecting_rod,
                            crank_mass1, connecting_rod_mass1, theta_values):
    """Shaking force after applying a balancing factor x (0 to 1) to the
    primary and secondary force components, simulating a counterweight
    that only partially cancels the reciprocating forces."""
    primary = ((piston_mass + connecting_rod_mass2) * crank * (w2 ** 2)) * np.cos(theta_values)
    secondary = ((piston_mass + connecting_rod_mass2) * (crank ** 2) * (w2 ** 2) / connecting_rod) * np.cos(2 * theta_values)
    fc = (crank_mass1 + connecting_rod_mass1) * crank * (w2 ** 2)
    fs = fc * np.exp(1j * theta_values) + x * primary + x * secondary
    return np.array([complex(force.evalf()) for force in fs], dtype=np.complex128)


def find_max_positive_negative(arr):
    """Find the largest positive value and the most negative value in an array."""
    if len(arr) == 0:
        return None, None

    max_positive = float('-inf')
    max_negative = float('inf')

    for num in arr:
        if num > 0 and num > max_positive:
            max_positive = num
        elif num < 0 and num < max_negative:
            max_negative = num

    if max_positive == float('-inf'):
        max_positive = None
    if max_negative == float('inf'):
        max_negative = None

    return max_positive, max_negative


# Search over candidate balancing factors x in [0, 1] and keep the ones
# that don't produce shaking forces larger than the unbalanced case
x_values = np.linspace(0, 1, 10)
in_range_x_values = []

for x in x_values:
    forces = forces_function(
        piston_mass, connecting_rod_mass2, crank, w2, connecting_rod,
        crank_mass1, connecting_rod_mass1, theta_values
    )
    max_positive_forces, max_negative_forces = find_max_positive_negative(forces)

    fs = shaking_force_function(
        x, piston_mass, connecting_rod_mass2, crank, w2, connecting_rod,
        crank_mass1, connecting_rod_mass1, theta_values
    )
    max_positive_fs, max_negative_fs = find_max_positive_negative(fs)

    if max_positive_fs < max_positive_forces and max_negative_fs < max_negative_forces:
        in_range_x_values.append(x)


def surface_function(force_values, theta_values):
    """Approximate the area enclosed by the force locus using the
    trapezoidal rule on the force magnitude vs. theta curve."""
    force_magnitudes = np.abs(force_values)
    delta_theta = np.diff(theta_values)
    surface_area = np.sum(0.5 * (force_magnitudes[:-1] + force_magnitudes[1:]) * delta_theta)
    return surface_area


# Among the acceptable balancing factors, find the one that minimizes
# the shaking-force surface area relative to the unbalanced case
final_x = []
safety_control = 0
minimum_surface = None

for x in in_range_x_values:
    forces = forces_function(
        piston_mass, connecting_rod_mass2, crank, w2, connecting_rod,
        crank_mass1, connecting_rod_mass1, theta_values
    )
    original_surface = surface_function(forces, theta_values)

    if safety_control == 0:
        minimum_surface = original_surface
        safety_control = 1

    fs = shaking_force_function(
        x, piston_mass, connecting_rod_mass2, crank, w2, connecting_rod,
        crank_mass1, connecting_rod_mass1, theta_values
    )
    fs_surface = surface_function(fs, theta_values)

    difference_surface = original_surface - fs_surface
    if 0 < difference_surface < minimum_surface:
        minimum_surface = difference_surface
        final_x = x

final_fs = shaking_force_function(
    final_x, piston_mass, connecting_rod_mass2, crank, w2, connecting_rod,
    crank_mass1, connecting_rod_mass1, theta_values
)

# Plot original vs. balanced (minimized) shaking forces
plt.figure(figsize=(10, 8))
plt.plot(np.real(force_values), np.imag(force_values), label='Original Forces', color='blue')
plt.plot(np.real(final_fs), np.imag(final_fs), label='Final Shaking Forces', linestyle='dashed', color='red')
plt.title('Original Forces and Final Shaking Forces')
plt.xlabel('Horizontal Forces')
plt.ylabel('Vertical Forces')
plt.legend()
plt.grid(True)
plt.show()

print(f"The suitable x value for our purpose is: {final_x}")
