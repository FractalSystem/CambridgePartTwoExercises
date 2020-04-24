import numpy as np
import scipy.integrate
import constants
import time


class Asteroid():
    def __init__(self, r_initial, v_initial, planet_mass = None):
        # Check if a custom planet mass is to be used for the system
        if planet_mass:
            self.MASS_JUPITER = planet_mass
        else:
            self.MASS_JUPITER = constants.MASS_JUPITER

        # Set initial conditions
        self.r_initial = r_initial
        self.v_initial = v_initial

        # Derive omega and orbital period
        self.omega = [0, 0, np.sqrt(constants.G * (constants.MASS_SUN + self.MASS_JUPITER) / constants.R ** 3)]
        self.orbital_period = 2 * np.pi / np.linalg.norm(self.omega)

        # Take coordinate centre as centre of mass of system
        self.r_s = np.array([-self.MASS_JUPITER * constants.R / (self.MASS_JUPITER + constants.MASS_SUN), 0, 0])  # Vector displacement from COM to Sun
        self.r_j = np.array([constants.MASS_SUN * constants.R / (self.MASS_JUPITER + constants.MASS_SUN), 0, 0])  # Vector displacement from COM to Jupiter

    def derivatives(self, t, y):
        r_a = y[:3]
        v_a = y[3:6]

        # Find displacement of asteroid from masses
        r_a_to_s = self.r_s - r_a
        r_a_to_j = self.r_j - r_a

        # Define force per unit asteroid mass due to gravity
        F = constants.G * constants.MASS_SUN * r_a_to_s / (np.linalg.norm(r_a_to_s) ** 3) + constants.G * self.MASS_JUPITER * r_a_to_j / (
                np.linalg.norm(r_a_to_j) ** 3)

        # Define equations of motion
        r_a_dot = v_a
        v_a_dot = F - 2 * np.cross(self.omega, v_a) - np.cross(self.omega, np.cross(self.omega, r_a))

        return np.concatenate((r_a_dot, v_a_dot))

    def derivatives2(self, t, y):
        r_x, r_y, r_z, v_x, v_y, v_z = y

        # Manually define displacements of asteroid from masses, without using Numpy arrays
        r_a_to_s_x = self.r_s[0] - r_x
        r_a_to_s_y = -r_y
        r_a_to_j_x = self.r_j[0] - r_x
        r_a_to_j_y = -r_y

        # Define moduli cubed to avoid double calculation
        mod_r_a_to_s_cubed = (r_a_to_s_x ** 2 + r_a_to_s_y ** 2) ** 1.5
        mod_r_a_to_j_cubed = (r_a_to_j_x ** 2 + r_a_to_j_y ** 2) ** 1.5

        # Find componenets of F seperately
        F_x = constants.G * constants.MASS_SUN * (r_a_to_s_x) / (
            mod_r_a_to_s_cubed) + constants.G * self.MASS_JUPITER * r_a_to_j_x / (
                  mod_r_a_to_j_cubed)
        F_y = constants.G * constants.MASS_SUN * (r_a_to_s_y) / (
            mod_r_a_to_s_cubed) + constants.G * self.MASS_JUPITER * r_a_to_j_y / (
                  mod_r_a_to_j_cubed)

        # Return as a tuple, without creating a new Numpy array
        return (
            v_x, v_y, v_z, F_x + 2 * self.omega[2] * v_y + r_x * self.omega[2] ** 2,
            F_y - 2 * self.omega[2] * v_x + r_y * self.omega[2] ** 2,
            0)

    def solve_orbit(self, n_orbits):
        # initial_conditions are given when class is initiated
        start_time = time.time()
        t_max = n_orbits * self.orbital_period
        solution = scipy.integrate.solve_ivp(self.derivatives2,
                                             t_span=(0, t_max),
                                             t_eval=np.linspace(0, t_max, 200*n_orbits),
                                             y0=np.concatenate((self.r_initial, self.v_initial)),
                                             rtol=10**-4,
                                             atol = 10**-7,
                                             method = scipy.integrate.RK45)

        t, r_a, v_a = (solution.t, solution.y[:3], solution.y[3:6])
        # print(f"Took {time.time()-start_time}s")
        return t, r_a, v_a
