#######################
# File containing the class Asteroid. Asteroid is used to represent a single asteroid's initial conditions and simulate
# its trajectory over a specified number of orbits.
#######################

import numpy as np
import scipy.integrate
import constants


class Asteroid:
    def __init__(self, r_initial, v_initial, planet_mass=None):
        """
        Constructor of Asteroid class. Defines initial conditions of asteroids
        :param r_initial: (np.ndarray) 3 dimensional array of initial position
        :param v_initial: (np.ndarray) 3 dimensional array of initial velocity
        :param planet_mass: (float) Optional value used to specify a custom Planet-Sun mass ratio
        """
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
        self.r_s = np.array([-self.MASS_JUPITER * constants.R / (self.MASS_JUPITER + constants.MASS_SUN), 0,
                             0])  # Vector displacement from COM to Sun
        self.r_j = np.array([constants.MASS_SUN * constants.R / (self.MASS_JUPITER + constants.MASS_SUN), 0,
                             0])  # Vector displacement from COM to Jupiter

    def derivatives(self, t, y):
        """
        Depreciated function using Numpy arrays and methods to calculate the coupled equations of motion. Slow compared
        to derivatives2. Called by scipy.integrate.solve_ivp().
        :param t: (np.float64) Scalar value of integration (time in this case).
        :param y: (np.ndarray) Array containing all input values at a step in the integration.
        :return: (np.ndarray) Array of same shape as y containing calculated results from the equations of motion.
        """
        r_a = y[:3]
        v_a = y[3:6]

        # Find displacement of asteroid from masses
        r_a_to_s = self.r_s - r_a
        r_a_to_j = self.r_j - r_a

        # Define force per unit asteroid mass due to gravity
        F = constants.G * constants.MASS_SUN * r_a_to_s / (
                    np.linalg.norm(r_a_to_s) ** 3) + constants.G * self.MASS_JUPITER * r_a_to_j / (
                    np.linalg.norm(r_a_to_j) ** 3)

        # Define equations of motion
        r_a_dot = v_a
        v_a_dot = F - 2 * np.cross(self.omega, v_a) - np.cross(self.omega, np.cross(self.omega, r_a))

        return np.concatenate((r_a_dot, v_a_dot))

    def derivatives2(self, t, y):
        """
        Function using element-wise operations to calculate the coupled equations of motion. Faster than equivalent
        Numpy method. Called by scipy.integrate.solve_ivp().
        :param t: (np.float64) Scalar value of integration (time in this case).
        :param y: (np.ndarray) Array containing all input values at a step in the integration.
        :return: (tuple) Tuple of same shape as y containing calculated results from the equations of motion.
        """
        # Unpack input variables
        r_x, r_y, r_z, v_x, v_y, v_z = y

        # Manually define displacements of asteroid from masses, without using Numpy arrays
        r_a_to_s_x = self.r_s[0] - r_x
        r_a_to_s_y = -r_y
        r_a_to_j_x = self.r_j[0] - r_x
        r_a_to_j_y = -r_y

        # Define moduli cubed to avoid double calculation
        mod_r_a_to_s_cubed = (r_a_to_s_x ** 2 + r_a_to_s_y ** 2) ** 1.5
        mod_r_a_to_j_cubed = (r_a_to_j_x ** 2 + r_a_to_j_y ** 2) ** 1.5

        # Find components of F separately
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
        """
        Function called to simulate trajectory of Trojan asteroid for n_orbits
        :param n_orbits: (int) Number of orbits of Jupiter to simulate the asteroid for.
        :return: ((np.ndarray), (np.ndarray), (np.ndarray)) Arrays of t, y and y' representing the solution of the
                    simulation.
        """
        # Define time in years up to which to simulate the asteroid path
        t_max = n_orbits * self.orbital_period

        # Use RK45 ODE solver to solve initial value problem. Adjust tolerances using rtol and atol.
        solution = scipy.integrate.solve_ivp(self.derivatives2,
                                             t_span=(0, t_max),
                                             t_eval=np.linspace(0, t_max, 200 * n_orbits),
                                             y0=np.concatenate((self.r_initial, self.v_initial)),
                                             rtol=10 ** -4,
                                             atol=10 ** -7,
                                             method=scipy.integrate.RK45)

        # Return the solution
        t, r_a, v_a = (solution.t, solution.y[:3], solution.y[3:6])
        return t, r_a, v_a
