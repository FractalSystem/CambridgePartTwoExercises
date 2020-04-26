#######################
# Small file which is loaded by each created process during evaluation of wander.
#######################


from time import time
from asteroid import Asteroid
import numpy as np


def pooled_process_position(args):
    """
    Evaluate wander based on input list of initial position values.
    :param args: (list) List containing initial conditions of the form (x, y, r_lagrange_point). Each item is to be
                    evaluated and have its wander calculated.
    :return: (list) List of calculated ranges of wander in the same order as the input list.
    """
    # Iterate through each item in args, calculate wander and append to r_max_list
    r_max_list = []
    for item in args:
        # Save start time of function
        time_start = time()

        # Unpack item and set initial conditions
        x, y, r_lagrange_point = item
        r_a_initial = [x, y, 0]
        v_a_initial = [0, 0, 0]

        # Simulate asteroid for 100 orbits
        asteroid = Asteroid(r_a_initial, v_a_initial)
        t, r_a, v_a = asteroid.solve_orbit(100)

        # Calculate r_max, the maximum distance of the asteroid from the Lagrange point
        r_max = np.amax(np.sqrt((r_a[0] - r_lagrange_point[0]) ** 2 + (r_a[1] - r_lagrange_point[1]) ** 2))
        r_max_list.append(r_max)

        # Output runtime for each loop to provide feedback during long calculation
        print(f"Took {time()-time_start}s")

    return r_max_list


def pooled_process_velocity(args):
    """
    Evaluate wander based on input list of initial velocity values.
    :param args: (list) List containing initial conditions of the form (v_x, v_y, r_lagrange_point). Each item is to be
                    evaluated and have its wander calculated.
    :return: (list) List of calculated ranges of wander in the same order as the input list.
    """
    # Iterate through each item in args, calculate wander and append to r_max_list
    r_max_list = []
    for item in args:
        # Save start time of function
        time_start = time()

        # Unpack item and set initial conditions
        x, y, r_lagrange_point = item
        r_a_initial = [r_lagrange_point[0], r_lagrange_point[1], 0]
        v_a_initial = [x, y, 0]

        # Simulate asteroid for 100 orbits
        asteroid = Asteroid(r_a_initial, v_a_initial)
        t, r_a, v_a = asteroid.solve_orbit(100)

        # Calculate r_max, the maximum distance of the asteroid from the Lagrange point
        r_max = np.amax(np.sqrt((r_a[0] - r_lagrange_point[0]) ** 2 + (r_a[1] - r_lagrange_point[1]) ** 2))
        r_max_list.append(r_max)

        # Output runtime for each loop to provide feedback during long calculation
        print(f"Took {time()-time_start}s")

    return r_max_list
