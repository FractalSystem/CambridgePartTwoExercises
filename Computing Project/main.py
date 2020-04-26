#######################
# File defines the class Main which contains all logic, plotting and saving functions. It also serves as the
# entry point of the program.
#######################

import numpy as np
import matplotlib.pyplot as plt
import constants
import matplotlib.animation
from multiprocessing import Pool
from asteroid import Asteroid
from pooled_process import pooled_process_position, pooled_process_velocity
import json
from time import time
import scipy.stats


class Main:
    def __init__(self):
        """
        Constructor method for class main. Used to set up vectors r_s and r_j to prevent multiple recalculations of them.
        """
        # Define vector displacement from COM to Sun
        self.r_s = np.array([-constants.MASS_JUPITER * constants.R / (constants.MASS_JUPITER + constants.MASS_SUN), 0,
                             0])

        # Define vector displacement from COM to Jupiter
        self.r_j = np.array([constants.MASS_SUN * constants.R / (constants.MASS_JUPITER + constants.MASS_SUN), 0,
                             0])

    def plot_extras(self, ax1):
        """
        Plots graphical extras on subplot defined by ax1. Extras are Jupiter, the Sun, COM and Jupiter's radius of orbit.
        :param ax1: (matplotlib.axes._subplots.AxesSubplot) Axis to have extras plotted on to.
        :return: None.
        """
        # Plot radius of orbit of jupiter
        orbit_circle = plt.Circle((self.r_s[0], self.r_s[1]), constants.R, fill=False, linewidth=0.5, linestyle='--')
        ax1.add_artist(orbit_circle)

        # Plot Sun's location
        sun_circle = plt.Circle((self.r_s[0], self.r_s[1]), 0.45, linewidth=0.3, color='#ffff4b', ec='k')
        ax1.add_artist(sun_circle)

        # Plot Jupiter's location
        jupiter_circle = plt.Circle((self.r_j[0], self.r_j[1]), 0.1, color="#9a9aff")
        ax1.add_artist(jupiter_circle)

        # Plot COM
        ax1.plot(0, 0, 'b+')

    def plot_orbit(self, r_a_initial, v_a_initial, n_orbits=100):
        """
        Plotting routine for a standard, two picture, orbital plot.
        :param r_a_initial: (list) Initial position vector of the asteroid.
        :param v_a_initial: (list) Initial velocity vector of the asteroid
        :param n_orbits: (int) Number of orbits of Jupiter to plot.
        :return: None
        """

        # Define asteroid and simulate its trajectory for n_orbits orbits of Jupiter
        asteroid = Asteroid(r_a_initial, v_a_initial)
        t, r_a, v_a = asteroid.solve_orbit(n_orbits)

        # Plot overview of orbit (axis 1)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[3.0 * 4, 2.8 * 2])
        ax1.set_xlim((-7, 7))
        ax1.set_ylim((-7, 7))
        ax1.plot(r_a[0], r_a[1])
        self.plot_extras(ax1)
        ax1.set_xlabel("x /AU")
        ax1.set_ylabel("y /AU")

        # Plot zoomed in view of asteroid orbit (axis 2)
        ax2.plot(r_a[0], r_a[1])
        ax2.set_xlabel("x /AU")
        ax2.set_ylabel("y /AU")

        # Plot starting position on axis 2
        ax2.plot(r_a[0][0], r_a[1][0], 'r+')

        # Save figure
        plt.savefig("fig.png")

    def animate(self):
        """
        Animation routine for a standard orbit, slightly perturbed from L4. Used to create "animation1.mp4".
        :return: None.
        """
        # Define initial conditions
        r_a_initial = [constants.R * (
                (constants.MASS_SUN - constants.MASS_JUPITER) / (constants.MASS_SUN + constants.MASS_JUPITER)) * np.cos(
            np.pi / 3), constants.R * np.sin(np.pi / 3), 0]
        v_a_initial = [0, 0, 0]
        r_a_initial = np.array(r_a_initial) + np.array([-0.001 * constants.R, +0.001 * constants.R, 0])

        # Create asteroid object and simulate for 60 orbits
        asteroid = Asteroid(r_a_initial, v_a_initial)
        t, r_a, v_a = asteroid.solve_orbit(60)

        # Create figure and subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[3.2 * 4, 2.8 * 2])

        # Set axis limits and label axes
        ax1.set_xlim((-7, 7))
        ax1.set_ylim((-7, 7))
        ax2.set_xlim((np.min(r_a[0]), np.max(r_a[0])))
        ax2.set_ylim((np.min(r_a[1]), np.max(r_a[1])))
        ax1.set_xlabel("x /AU")
        ax1.set_ylabel("y /AU")
        ax2.set_xlabel("x /AU")
        ax2.set_ylabel("y /AU")

        # Plot Sun and Jupiter
        self.plot_extras(ax1)

        # Define variables to plot asteroid loci
        line, = ax1.plot([], [])
        line_zoomed, = ax2.plot([], [])
        point, = ax1.plot([], [], "ro", markersize=3)

        # Define graphic to show orbit location
        sun_circle = plt.Circle((-6, 6), 0.15, color="k")
        radius_line = plt.Circle((-6, 6), 0.8, color="k", fill=False, linewidth=0.2)
        ax1.add_artist(sun_circle)
        ax1.add_artist(radius_line)
        text = ax1.text(-4.8, 6, "t = 0 Yr")

        # Define orbit properties for animation
        omega = np.sqrt(constants.G * (constants.MASS_SUN + constants.MASS_JUPITER) / constants.R ** 3)
        period = 2 * np.pi / omega
        initial_angle = np.pi / 2  # Additional phase to make graphic plot match large plot at t=0

        # Define animation function
        def animate(i):
            line.set_data(r_a[0][0:i], r_a[1][0:i])
            line_zoomed.set_data(r_a[0][0:i], r_a[1][0:i])
            time = t[i]
            angle = ((time % period) / period) * 2 * np.pi + initial_angle
            x = -6 + 0.8 * np.sin(angle)
            y = 6 + 0.8 * np.cos(angle)
            point.set_data(x, y)
            text.set_text(f"t = {int(t[i])} Yr")
            return line, line_zoomed, point, text

        # Calculate interval and number of frames needed
        FPS = 60.0  # Frames per second in Hz
        ANIM_LENGTH = 20.0  # Animation length in seconds
        interval = 1 / FPS
        frames = int(ANIM_LENGTH * FPS)

        # Define animation
        animation = matplotlib.animation.FuncAnimation(fig, animate, frames=frames, interval=interval, blit=True)

        # Save animation
        plt.rcParams['animation.ffmpeg_path'] = constants.FFMPEG_PATH
        FFWriter = matplotlib.animation.writers['ffmpeg']
        writer = FFWriter(fps=FPS, metadata=dict(artist='Cambridge Computing Project 2020'), bitrate=2000)
        animation.save('animation1.mp4', writer=writer)

    def evaluate_wander(self, grid_size):
        """
        Function to evaluate wander over a grid of initial positions, and save the results.
        :param grid_size: (int) Used to define the dimensions of the grid over which to find wander.
        :return: None.
        """

        #

        # Define range of values in each direction to add to L4 position
        r_values = np.linspace(-2, +2, grid_size)

        # Define position of L4
        r_lagrange_point = [constants.R * (
                (constants.MASS_SUN - constants.MASS_JUPITER) / (constants.MASS_SUN + constants.MASS_JUPITER)) * np.cos(
            np.pi / 3), constants.R * np.sin(np.pi / 3), 0]

        # Define 2D arrays X, Y about the lagrange point for the wander to be evaluated at
        X, Y = np.meshgrid(r_values + r_lagrange_point[0], r_values + r_lagrange_point[1])

        # Generate input_list to supply initial conditions to worker processes
        input_list = []
        for i in range(grid_size):
            for j in range(grid_size):
                input_list.append([X[i][j], Y[i][j], r_lagrange_point])

        # Split input list into sections of length n
        n = int(len(input_list) / 4)
        input_list = [input_list[i:i + n] for i in range(0, len(input_list), n)]

        # Define pool and map input_list to the pooled processes
        pool = Pool()
        result = pool.map(pooled_process_position, input_list)

        # 'result' is a 2D array (i, j) corresponding to coordinates X[i][j], Y[i][j]
        result = np.concatenate(result).reshape((grid_size, grid_size))

        # Save result to disk
        self.save_results(result, X, Y)

    def evaluate_wander_velocity(self, grid_size):
        """
        Function to evaluate wander over a grid of initial velocities, and save the results.
        :param grid_size: (int) Used to define the dimensions of the grid over which to find wander.
        :return: None.
        """

        # Define range of velocities either side of zero
        v_values = np.linspace(-0.4, +0.4, grid_size)

        # Define position of L4
        r_lagrange_point = np.array([constants.R * np.sin(np.pi / 6),
                                     constants.R * ((constants.MASS_SUN - constants.MASS_JUPITER) / (
                                             constants.MASS_SUN + constants.MASS_JUPITER)) * np.cos(np.pi / 6),
                                     0])  # Asteroid vector displacement from COM

        # Define 2D arrays X, Y about the lagrange point for the wander to be evaluated at
        X, Y = np.meshgrid(v_values, v_values)

        # Generate input_list to supply initial conditions to worker processes
        input_list = []
        for i in range(grid_size):
            for j in range(grid_size):
                input_list.append([X[i][j], Y[i][j], r_lagrange_point])

        # Split input list into sections of length n
        n = int(len(input_list) / 4)
        input_list = [input_list[i:i + n] for i in range(0, len(input_list), n)]

        # Define pool and map input_list to the pooled processes
        pool = Pool()
        result = pool.map(pooled_process_velocity, input_list)

        # 'result' is a 2D array (i, j) corresponding to coordinates X[i][j], Y[i][j]
        result = np.concatenate(result).reshape((grid_size, grid_size))

        # Save result to disk
        self.save_results(result, X, Y)

    def save_results(self, result, X, Y):
        """
        Small function to save results from wander evaluations to disk.
        :param result: (np.ndarray) result from wander evaluation.
        :param X: (np.ndarray) X array from meshgrid().
        :param Y: (np.ndarray) Y array from meshgrid().
        :return: None.
        """
        # Define dictionary to be serialized to JSON
        dic = {"X": X.tolist(), "Y": Y.tolist(), "results": result.tolist()}

        # Open file in write mode and dump dictionary to JSON string
        with open("results.txt", 'w') as f:
            f.write(json.dumps(dic))

    def load_results(self, filename):
        """
        Small function to load results from wander evaluations from disk.
        :param filename: (str) Filename to load.
        :return: ((np.ndarray), (np.ndarray), (np.ndarray)) Returns X, Y, result, as defined in save_results.
        """
        # Open file in read mode and load JSON strong to dictionary object
        with open(filename, "r") as f:
            dic = json.loads(f.read())
        return dic.get("X"), dic.get("Y"), dic.get("results")

    def plot_wander(self, X, Y, results):
        """
        Plotting function for wander evaluation. Usually called after loading results from disk.
        :param X: (np.ndarray) X array from meshgrid().
        :param Y: (np.ndarray) Y array from meshgrid().
        :param result: (np.ndarray) result from wander evaluation.
        :return: None.
        """
        # Define figure
        fig, ax1 = plt.subplots()

        # Define position of L4 and plot as blue cross.
        r_a_initial = [constants.R * (
                (constants.MASS_SUN - constants.MASS_JUPITER) / (constants.MASS_SUN + constants.MASS_JUPITER)) * np.cos(
            np.pi / 3), constants.R * np.sin(np.pi / 3), 0]
        ax1.plot(r_a_initial[0], r_a_initial[1], "b+")

        # Plot contour map of log(results) with 1000 contours
        cp = ax1.contourf(X, Y, np.log(results) / np.log(10), 1000)

        # Label axes and colour bar scale
        ax1.set_xlabel("$v_x$ /AUYr$^{-1}$")
        ax1.set_ylabel("$v_y$ /AUYr$^{-1}$")
        cbar = fig.colorbar(cp)
        cbar.set_label('log$_{10}$(Wander /AU)')

        # Save figure
        plt.savefig("fig.png")

    def plot_potential(self, mass_ratio):
        """
        Evaluate and plot potential energy contours by evaluating over a grid of coordinates.
        :param mass_ratio: (float) Decimal ratio defined by MASS_PLANET/MASS_SUN. Used to customise planet mass for clearer plot.
        :return: None.
        """
        # Modify ratio of masses to get illustrative plot
        MODIFIED_MASS_JUPITER = mass_ratio * constants.MASS_SUN

        # Re-derive r_s for modified mass ratio
        r_s = np.array([-MODIFIED_MASS_JUPITER * constants.R / (MODIFIED_MASS_JUPITER + constants.MASS_SUN), 0,
                        0])

        # Derive angular velocity of rotating frame
        omega = np.sqrt(constants.G * (constants.MASS_SUN + MODIFIED_MASS_JUPITER) / constants.R ** 3)

        # Initiate 2D grid centred on the sun
        xlist = np.linspace(-7.0, 7.0, 1000)
        ylist = np.linspace(-7.0, 7.0, 1000)
        X, Y = np.meshgrid(xlist, ylist)

        # Effective potential is sum of potential from each mass + centripetal force potential. Coriolis force cannot be plotted due to dependence on asteroid velocity.
        potential = -constants.G * (constants.MASS_SUN / np.sqrt(X ** 2 + Y ** 2) + MODIFIED_MASS_JUPITER / np.sqrt(
            (X - constants.R) ** 2 + Y ** 2)) - 1 / 2 * omega ** 2 * ((X + r_s[0]) ** 2 + Y ** 2)

        # Find and mark positions of maxima in potential
        Z_max = np.amax(potential)
        maxima = np.where(np.isclose(potential, Z_max, rtol=10 ** -8))  # maxima is tuple : ([x_coords], [y_coords])
        fig, ax1 = plt.subplots()
        for i in range(len(maxima[0])):
            x_coord, y_coord = maxima[0][i], maxima[1][i]
            x, y = X[x_coord][y_coord], Y[x_coord][y_coord]
            ax1.plot(x, y, 'b+')

        # Plot contour plot of potential
        cp = ax1.contour(X, Y, potential, np.linspace(Z_max - 8, Z_max - 0.004, 200))

        # Label colour bar and axes
        cbar = fig.colorbar(cp)
        cbar.set_label('Potential /M$_{\odot}$AU$^2$Yr$^{-2}$')
        ax1.set_xlabel("x /AU")
        ax1.set_ylabel("y /AU")

        # Save figure
        plt.savefig("fig.png")

    def evaluate_mass_wander(self):
        """
        Evaluate wander of an asteroid from L4 for a range of planetary masses.
        :return: None.
        """
        # Define initial velocity vector
        v_a_initial = [0, 0, 0]

        # Populate an array of masses to iterate over
        # masses = np.linspace(0.0001, 0.05, num=100) * constants.MASS_SUN
        masses = np.logspace(-5, -2, num=100)

        # Evaluate masses between a given range of solar masses.
        r_max_list = []
        for mass in masses:
            # Recalculate Lagrange point (L4) position
            r_lagrange = [constants.R * ((constants.MASS_SUN - mass) / (
                    constants.MASS_SUN + mass)) * np.cos(np.pi / 3), constants.R * np.sin(np.pi / 3), 0]

            # Perturb Lagrange point
            # r_lagrange = np.array(r_lagrange) + np.array([-0.0006 * constants.R, +0.0006 * constants.R, 0])

            # Define a new asteroid with specified planet mass
            asteroid = Asteroid(r_lagrange, v_a_initial, planet_mass=mass)

            # Simulate asteroid over 500 orbits of Jupiter
            t, r_a, v_a = asteroid.solve_orbit(500)

            # Find range of wander r_max
            r_max = np.amax(np.sqrt((r_a[0] - r_lagrange[0]) ** 2 + (r_a[1] - r_lagrange[1]) ** 2))

            # Append range of wander to list to be plotted
            r_max_list.append(r_max)

            # Print output as loop executes for debugging
            print(r_max, mass)

        # Open file and dump results to it for future evaluation
        with open("out.txt", "w") as f:
            f.write(json.dumps({"r_max_list": r_max_list, "masses": masses.tolist()}))

    def plot_mass_wander(self, max_fit_index=None):
        """
        Plot and perform linear regression on previously calculated mass ratio - wander data, stored in the file "out.txt".
        :param max_fit_index: (int) Index up to which to include in linear regression analysis.
        :return: None
        """

        # Load data from disk
        with open("out.txt", "r") as f:
            dic = json.loads(f.read())
        r_max_list = dic.get("r_max_list")
        masses = dic.get("masses")

        # If no max_fit_index then just perform linear regression on entire data-set
        if max_fit_index is None:
            max_fit_index = len(masses)

        # Take natural log of masses and r_max_list for second graph and fitting purposes
        ln_masses = np.log(masses)
        ln_r_max = np.log(r_max_list)

        # Fit log-log graph using linear regression
        slope, intercept, r_val, p_val, std_err = scipy.stats.linregress(ln_masses[:max_fit_index],
                                                                         ln_r_max[:max_fit_index])
        fit_coefficients = [slope, intercept]

        # Output linear fit results
        print(
            f"Linear regression of log-log found relationship: y = kx^n where k = {np.exp(fit_coefficients[1])} and n = {fit_coefficients[0]}, with an R^2 value of {r_val**2}")

        # Define plotting parameters
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[3.0 * 4, 2.8 * 2])

        # Plot data
        ax1.plot(masses, r_max_list, label="Data")
        ax1.plot(masses, np.exp(fit_coefficients[1]) * masses ** fit_coefficients[0], label="Fit")
        ax2.plot(ln_masses, ln_r_max, label="Data", marker="+")
        ax2.plot(ln_masses, fit_coefficients[0] * ln_masses + fit_coefficients[1], label="Fit")

        # Set axis labels
        ax1.set_xlabel("Planet-Sun mass ratio /No Units")
        ax1.set_ylabel("Range of wander /AU")
        ax2.set_xlabel("ln(Planet-Sun mass ratio /No Units)")
        ax2.set_ylabel("ln(Range of wander /AU)")

        # Include legends
        ax1.legend()
        ax2.legend()

        # Save figure
        plt.savefig("fig.png")

    def evaluate_energy_conservation(self, r_a_initial, v_a_initial, n_orbits):
        """
        Test energy conservation for a given number of orbits and initial conditions.
        :param r_a_initial: (list) Initial position vector of the asteroid.
        :param v_a_initial: (list) Initial velocity vector of the asteroid
        :param n_orbits: (int) Number of orbits of Jupiter to plot.
        :return: None.
        """
        # Define asteroid object with initial conditions
        asteroid = Asteroid(r_a_initial, v_a_initial)

        # Simulate orbit for n_orbits
        t, r_a, v_a = asteroid.solve_orbit(n_orbits)

        # Transpose r_a and v_a for future convenience
        r_a = np.transpose(r_a)
        v_a = np.transpose(v_a)

        # Define omega and orbital period
        omega = np.sqrt(constants.G * (constants.MASS_SUN + constants.MASS_JUPITER) / constants.R ** 3)
        period = 2 * np.pi / omega

        # Iterate over each time sample and calculate total energy
        energy_list = []
        for i in range(len(t)):
            # Define vectors for potential calculation
            r_a_to_s = self.r_s - r_a[i]
            r_a_to_j = self.r_j - r_a[i]
            r = np.linalg.norm(r_a[i])

            # Calculate gravitational potential energy
            graviatational_potential = -constants.G * (
                    constants.MASS_SUN / np.linalg.norm(r_a_to_s) + constants.MASS_JUPITER / np.linalg.norm(r_a_to_j))

            # Calculate kinetic energy due to rotating frame
            kinetic_energy_rot = 0.5 * r ** 2 * omega ** 2

            # Calculate kinetic energy in the rotating frame
            kinetic_energy_frame = 0.5 * np.linalg.norm(v_a[i]) ** 2

            # Calculate total kinetic energy
            kinetic_energy = -kinetic_energy_rot + kinetic_energy_frame

            # Append total energy to energy_list
            energy_list.append(kinetic_energy + graviatational_potential)

        fig, ax1 = plt.subplots()

        # Plot all energies as a percentage deviation from initial energy
        ax1.plot(t / period, (np.array(energy_list) / energy_list[0] - 1) * 100)

        # Plot expected theoretical energy line
        ax1.hlines(0, np.min(t), np.max(t / period))

        # Set axis labels
        ax1.set_xlabel("Number of orbits /No Units")
        ax1.set_ylabel("|Percentage deviation from expected energy| /%")

        # Save plot
        plt.savefig("fig.png")

    def time_it(self, func, args):
        """
        Small function to time another function call.
        :param func: (method) Function name to call.
        :param args: (tuple) Arguments to pass into function call.
        :return: None.
        """
        # Record start_time
        start_time = time()

        # Call the passed function, expanding args using splat operator *
        func(*args)

        # Print runtime
        runtime = time() - start_time
        print(f"Function call took {runtime}s")


if __name__ == "__main__":
    # Define main_obj used to call all functions in Main
    main_obj = Main()

    # Set initial conditions of asteroid
    r_a_initial = [constants.R * (
            (constants.MASS_SUN - constants.MASS_JUPITER) / (constants.MASS_SUN + constants.MASS_JUPITER)) * np.cos(
        np.pi / 3), constants.R * np.sin(np.pi / 3), 0]

    v_a_initial = [0, 0, 0]

    # Plot orbit from initial conditions
    main_obj.plot_orbit(r_a_initial, v_a_initial, n_orbits=285)

    # Create animation
    main_obj.animate()

    # Plot potential contour plot
    main_obj.plot_potential(0.01)

    # Evaluate wander using position phase space grid of 64x64
    main_obj.evaluate_wander(64)

    # Evaluate wander using velocity phase space grid of 64x64
    main_obj.evaluate_wander_velocity(64)

    # Load results from wander calculation and plot
    X, Y, results = main_obj.load_results("results.txt")
    main_obj.plot_wander(X, Y, results)

    # Evaluate wander as a function of Planet-Sun mass and plot from saved results
    main_obj.evaluate_mass_wander()
    main_obj.plot_mass_wander()
