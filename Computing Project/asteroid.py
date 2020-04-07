import numpy as np
import scipy.integrate
import constants



class Asteroid():
    def __init__(self, r_initial, v_initial):
        self.r_initial = r_initial
        self.v_initial = v_initial
        self.omega = [0, 0, np.sqrt(constants.G * (constants.MASS_SUN + constants.MASS_JUPITER) / constants.R ** 3)]
        self.orbital_period = 2 * np.pi / np.linalg.norm(self.omega)

        # Take coordinate centre as centre of mass of system
        self.r_s = np.array([-constants.MASS_JUPITER * constants.R / (constants.MASS_JUPITER + constants.MASS_SUN), 0, 0])  # Vector displacement from COM to Sun
        self.r_j = np.array([constants.MASS_SUN * constants.R / (constants.MASS_JUPITER + constants.MASS_SUN), 0, 0])  # Vector displacement from COM to Jupiter

    def derivatives(self, t, y):
        r_a = y[:3]
        v_a = y[3:6]

        # Find displacement of asteroid from masses
        r_a_to_s = self.r_s - r_a
        r_a_to_j = self.r_j - r_a

        # Define force per unit asteroid mass due to gravity
        F = constants.G * constants.MASS_SUN * r_a_to_s / (np.linalg.norm(r_a_to_s) ** 3) + constants.G * constants.MASS_JUPITER * r_a_to_j / (
                np.linalg.norm(r_a_to_j) ** 3)

        # Define equations of motion
        r_a_dot = v_a
        v_a_dot = F - 2 * np.cross(self.omega, v_a) - np.cross(self.omega, np.cross(self.omega, r_a))

        return np.concatenate((r_a_dot, v_a_dot))

    def solve_orbit(self, n_orbits):
        # initial_conditions are given when class is initiated
        t_max = n_orbits * self.orbital_period
        solution = scipy.integrate.solve_ivp(self.derivatives,
                                             t_span=(0, t_max),
                                             t_eval=np.linspace(0, t_max, 5000),
                                             y0=np.concatenate((self.r_initial, self.v_initial)))
        t, r_a, v_a = (solution.t, solution.y[:3], solution.y[3:6])
        return t, r_a, v_a