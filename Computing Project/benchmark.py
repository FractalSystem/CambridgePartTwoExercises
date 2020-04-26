#######################
# File responsible for benchmarking of derivatives functions
#######################

import numpy as np
import constants
import time

# Define constants
omega = [0, 0, np.sqrt(constants.G * (constants.MASS_SUN + constants.MASS_JUPITER) / constants.R ** 3)]
orbital_period = 2 * np.pi / np.linalg.norm(omega)
r_s = np.array([-constants.MASS_JUPITER * constants.R / (constants.MASS_JUPITER + constants.MASS_SUN), 0,
                0])
r_j = np.array([constants.MASS_SUN * constants.R / (constants.MASS_JUPITER + constants.MASS_SUN), 0,
                0])


def derivatives(t, y):
    # First function to benchmark. Uses Numpy arrays.
    r_a = y[:3]
    v_a = y[3:6]

    # Find displacement of asteroid from masses
    r_a_to_s = r_s - r_a
    r_a_to_j = r_j - r_a

    # Define force per unit asteroid mass due to gravity
    F = constants.G * constants.MASS_SUN * r_a_to_s / (
            np.linalg.norm(r_a_to_s) ** 3) + constants.G * constants.MASS_JUPITER * r_a_to_j / (
                np.linalg.norm(r_a_to_j) ** 3)

    # Define equations of motion
    r_a_dot = v_a
    v_a_dot = F - 2 * np.cross(omega, v_a) - np.cross(omega, np.cross(omega, r_a))

    return np.concatenate((r_a_dot, v_a_dot))


def derivatives2(t, y):
    # Second function to benchmark. Evaluates element-wise.
    r_x, r_y, r_z, v_x, v_y, v_z = y

    # Manually define displacements of asteroid from masses, without using Numpy arrays
    r_a_to_s_x = r_s[0] - r_x
    r_a_to_s_y = -r_y
    r_a_to_j_x = r_j[0] - r_x
    r_a_to_j_y = -r_y

    # Define moduli cubed to avoid double calculation
    mod_r_a_to_s_cubed = (r_a_to_s_x ** 2 + r_a_to_s_y ** 2) ** 1.5
    mod_r_a_to_j_cubed = (r_a_to_j_x ** 2 + r_a_to_j_y ** 2) ** 1.5

    # Find componenets of F seperately
    F_x = constants.G * constants.MASS_SUN * (r_a_to_s_x) / (
        mod_r_a_to_s_cubed) + constants.G * constants.MASS_JUPITER * r_a_to_j_x / (
              mod_r_a_to_j_cubed)
    F_y = constants.G * constants.MASS_SUN * (r_a_to_s_y) / (
        mod_r_a_to_s_cubed) + constants.G * constants.MASS_JUPITER * r_a_to_j_y / (
              mod_r_a_to_j_cubed)

    # Return as a tuple, without creating a new Numpy array
    return (
        v_x, v_y, v_z, F_x + 2 * omega[2] * v_y + r_x * omega[2] ** 2, F_y - 2 * omega[2] * v_x + r_y * omega[2] ** 2,
        0)


if __name__ == '__main__':
    # Define number of runs to benchmark over
    n_runs = 100000

    # For each run define some random starting conditions and time each function
    time_list_1 = []
    time_list_2 = []
    time_list_3 = []
    for n in range(n_runs):
        r_a = np.random.rand(2)
        v_a = np.random.rand(2)
        y1 = np.concatenate((r_a, [0], v_a, [0]))
        y2 = np.concatenate((r_a[0:2], v_a[0:2]))
        start_time = time.time()
        r = derivatives(1, y1)
        if n_runs == 1:
            print(r)
        run_time = time.time() - start_time
        time_list_1.append(run_time)
        start_time = time.time()
        r = derivatives2(1, y1)
        if n_runs == 1:
            print(r)
        run_time = time.time() - start_time
        time_list_2.append(run_time)

    # Print average runtimes
    print(np.average(time_list_1))
    print(np.average(time_list_2))
