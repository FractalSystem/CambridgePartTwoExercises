# Small file which is loaded by each created process during evaluation of wander.
from time import time
from asteroid import Asteroid
import numpy as np

def pooled_process_position(args):
    """
    Evaluate wander based on input list of initial values.
    :param args:
    :return:
    """
    r_max_list = []
    for item in args:
        time_start = time()
        x, y, r_lagrange_point = item
        r_a_initial = [x, y, 0]
        v_a_initial = [0, 0, 0]
        asteroid = Asteroid(r_a_initial, v_a_initial)
        t, r_a, v_a = asteroid.solve_orbit(100)

        # Calculate r_max, the maximum distance of the asteroid from the Lagrange point
        r_max = np.amax(np.sqrt((r_a[0] - r_lagrange_point[0]) ** 2 + (r_a[1] - r_lagrange_point[1]) ** 2))
        print(f"Took {time()-time_start}s")
        r_max_list.append(r_max)
    return r_max_list

def pooled_process_velocity(args):
    r_max_list = []
    for item in args:
        time_start = time()
        x, y, r_lagrange_point = item
        r_a_initial = [r_lagrange_point[0], r_lagrange_point[1], 0]
        v_a_initial = [x, y, 0]
        asteroid = Asteroid(r_a_initial, v_a_initial)
        t, r_a, v_a = asteroid.solve_orbit(100)

        # Calculate r_max, the maximum distance of the asteroid from the Lagrange point
        r_max = np.amax(np.sqrt((r_a[0] - r_lagrange_point[0]) ** 2 + (r_a[1] - r_lagrange_point[1]) ** 2))
        print(f"Took {time()-time_start}s")
        r_max_list.append(r_max)
    return r_max_list
