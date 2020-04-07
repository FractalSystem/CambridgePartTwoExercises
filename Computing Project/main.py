import numpy as np
import matplotlib.pyplot as plt
import constants
import matplotlib.animation
from multiprocessing import Pool
from asteroid import Asteroid
from pooled_process import pooled_process


### TEST ONLY
# r_a = [2.5948051948051956, 4.50333209967908, 0]
# omega = [0, 0,0.5298767365821111]



class Main():
    def __init__(self):
        # Take coordinate centre as centre of mass of system
        self.r_s = np.array([-constants.MASS_JUPITER * constants.R / (constants.MASS_JUPITER + constants.MASS_SUN), 0,
                             0])  # Vector displacement from COM to Sun

        self.r_j = np.array([constants.MASS_SUN * constants.R / (constants.MASS_JUPITER + constants.MASS_SUN), 0,
                             0])  # Vector displacement from COM to Jupiter

    def plot_extras(self, ax1):
        # Plot radius of orbit of jupiter
        orbit_circle = plt.Circle((self.r_s[0], self.r_s[1]), constants.R, fill=False, linewidth=0.5, linestyle='--')
        ax1.add_artist(orbit_circle)

        # Plot Sun's location
        sun_circle = plt.Circle((self.r_s[0], self.r_s[1]), 0.45, linewidth=0.3, color='#ffff4b', ec='k')
        ax1.add_artist(sun_circle)

        # Plot Jupiter's location
        jupiter_circle = plt.Circle((self.r_j[0], self.r_j[1]), 0.1, color = "#9a9aff")
        ax1.add_artist(jupiter_circle)

        # Plot COM
        ax1.plot(0, 0, 'b+')

    def plot_orbit(self):
        # Gives the default orbital view
        # Initial conditions of asteroid
        r_a_initial = [constants.R * np.sin(np.pi / 6),
                       constants.R * ((constants.MASS_SUN - constants.MASS_JUPITER) / (
                                   constants.MASS_SUN + constants.MASS_JUPITER)) * np.cos(
                           np.pi / 6),
                       0]  # Asteroid vector displacement from COM
        r_a_initial = np.array(r_a_initial) + np.array([-0.05, +0.05, 0])  # CARE! Perturbing initial radius
        v_a_initial = [0, 0, 0]

        # Define asteroid and simulate its orbit
        asteroid = Asteroid(r_a_initial, v_a_initial)
        t, r_a, v_a = asteroid.solve_orbit(100)

        # Plot overview of orbit
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[3.2 * 4, 2.8 * 2])
        ax1.set_xlim((0, 7))
        ax1.set_ylim((0, 7))
        ax1.plot(r_a[0], r_a[1])
        self.plot_extras(ax1)
        ax1.set_xlabel("x /Au")
        ax1.set_ylabel("y /Au")

        # Plot zoomed in view of asteroid orbit
        ax2.plot(r_a[0], r_a[1])
        ax2.set_xlabel("x /Au")
        ax2.set_ylabel("y /Au")

        # plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    def animate(self):

        # Model orbit of asteroid about Lagrange point
        r_a_initial = [constants.R * np.sin(np.pi / 6),
                       constants.R * ((constants.MASS_SUN - constants.MASS_JUPITER) / (
                                   constants.MASS_SUN + constants.MASS_JUPITER)) * np.cos(
                           np.pi / 6),
                       0]  # Asteroid vector displacement from COM
        v_a_initial = [0, 0, 0]
        asteroid = Asteroid(r_a_initial, v_a_initial)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[3.2 * 4, 2.8 * 2])
        t, r_a, v_a = asteroid.solve_orbit(60)

        # Set axes limits and label axes
        ax1.set_xlim((-7, 7))
        ax1.set_ylim((-7, 7))
        ax2.set_xlim((np.min(r_a[0]), np.max(r_a[0])))
        ax2.set_ylim((np.min(r_a[1]), np.max(r_a[1])))
        ax1.set_xlabel("x /Au")
        ax1.set_ylabel("y /Au")
        ax2.set_xlabel("x /Au")
        ax2.set_ylabel("y /Au")

        # Plot Sun and Jupiter
        self.plot_extras(ax1)

        # Plot asteroid loci
        line, = ax1.plot([], [])
        line_zoomed, = ax2.plot([], [])
        point, = ax1.plot([], [], "ro", markersize=3)

        # Define graphic to show orbit location
        sun_circle = plt.Circle((6, -6), 0.15, color="k")
        radius_line = plt.Circle((6, -6), 0.8, color="k", fill=False, linewidth=0.2)
        ax1.add_artist(sun_circle)
        ax1.add_artist(radius_line)

        # Define orbit properties for graphic
        omega = np.sqrt(constants.G * (constants.MASS_SUN + constants.MASS_JUPITER) / constants.R ** 3)
        period = 2*np.pi/omega
        initial_angle = np.pi/2 # Additional phase to make graphic plot match large plot at t=0


        # Define animation function
        def animate(i):
            line.set_data(r_a[0][0:i], r_a[1][0:i])
            line_zoomed.set_data(r_a[0][0:i], r_a[1][0:i])
            time = t[i]
            angle = ((time%period)/period)*2*np.pi+initial_angle
            x = 6+0.8*np.sin(angle)
            y = -6+0.8*np.cos(angle)
            point.set_data(x,y)
            return line, line_zoomed, point

        # Calculate interval and number of frames needed
        FPS = 60.0 # Frames per second in Hz
        ANIM_LENGTH = 20.0  # Animation length in seconds
        interval = 1/FPS
        frames = int(ANIM_LENGTH*FPS)

        # Define animation
        animation = matplotlib.animation.FuncAnimation(fig, animate, frames=frames, interval=interval, blit=True)

        # Save animation
        plt.rcParams['animation.ffmpeg_path'] = constants.FFMPEG_PATH
        FFWriter = matplotlib.animation.writers['ffmpeg']
        writer = FFWriter(fps=FPS, metadata=dict(artist='Noah Crew-Gee (nc506)'), bitrate=2000)
        animation.save('out.mp4', writer=writer)

        # plt.show()

    def evaluate_wander(self):
        # Only for initial position for now
        # Define a set of initial conditions in r_a, v_a phase space
        # Call solve orbit in parallel
        # Evaluate maximum wander using R_max
        grid_size = 4
        r_values = np.linspace(-0.05, +0.05, grid_size)
        r_lagrange_point = np.array([constants.R * np.sin(np.pi / 6),
                                     constants.R * ((constants.MASS_SUN - constants.MASS_JUPITER) / (
                                             constants.MASS_SUN + constants.MASS_JUPITER)) * np.cos(np.pi / 6),
                                     0])  # Asteroid vector displacement from COM

        r_lagrange_point = [5,20,0]
        X, Y = np.meshgrid(r_values+r_lagrange_point[0], r_values+r_lagrange_point[1])

        print(X)
        print(Y)
        from time import sleep
        sleep(3)

        r_max_array = np.zeros((grid_size, grid_size))

        # Generate input_list to supply initial conditions to worker processes
        input_list = []
        for i in range(grid_size):
            for j in range(grid_size):
                input_list.append((X, Y, r_lagrange_point, i, j))

        # Split input list into sections of length n
        n = int(len(input_list)/4)
        input_list = [input_list[i:i + n] for i in range(0, len(input_list), n)]

        # Define pool and map input_list to the pooled processes
        pool = Pool(4)
        result = pool.map(pooled_process, input_list)
        print(result)
        print(np.concatenate(result))


    def plot_potential(self):
        # Derive angular velocity of rotating frame
        omega = np.sqrt(constants.G * (constants.MASS_SUN + constants.MASS_JUPITER) / constants.R ** 3)

        # Initiate 2D grid centred on the sun
        xlist = np.linspace(-7.0, 7.0, 1000)
        ylist = np.linspace(-7.0, 7.0, 1000)
        X, Y = np.meshgrid(xlist, ylist)

        # Effective potential is sum of potential from each mass + centripetal force potential. Coriolis force cannot be plotted due to dependence on asteroid velocity.
        potential = -constants.G * (constants.MASS_SUN / np.sqrt(X ** 2 + Y ** 2) + constants.MASS_JUPITER / np.sqrt(
            (X - constants.R) ** 2 + Y ** 2)) - 1 / 2 * omega ** 2 * ((X + self.r_s[0]) ** 2 + Y ** 2)

        # potential from com has same problem
        potential = -constants.G * (constants.MASS_SUN/ np.sqrt((X-self.r_s[0])**2 + Y**2)+constants.MASS_JUPITER/
                                    np.sqrt((X-self.r_j[0])**2+Y**2)) - 1/2 * omega**2 * (X**2+Y**2)


        # Find and mark positions of maxima in Z
        Z_max = np.amax(potential)
        maxima = np.where(np.isclose(potential, Z_max, rtol=10 ** -8))  # maxima is tuple : ([x_coords], [y_coords])
        fig, ax1 = plt.subplots(1, 1)
        for i in range(len(maxima[0])):
            x_coord, y_coord = maxima[0][i], maxima[1][i]
            x, y = X[x_coord][y_coord], Y[x_coord][y_coord]
            ax1.plot(x, y, 'b+')

        # Plot contour plot of potential
        cp = ax1.contour(X, Y, potential, np.linspace(Z_max - 8, Z_max - 0.004, 200))
        fig.colorbar(cp)
        ax1.set_xlabel("x /Au")
        ax1.set_ylabel("y /Au")
        plt.show()


if __name__ == "__main__":
    main_obj = Main()
    # main_obj.plot_orbit()
    # main_obj.animate()
    main_obj.evaluate_wander()
    # main_obj.plot_potential()
