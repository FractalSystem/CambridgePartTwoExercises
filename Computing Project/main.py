import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
from multiprocessing import Pool

# Define constants
G = 4 * np.pi ** 2
mass_jupiter = 0.001  # in solar masses
mass_sun = 1.0  # in solar masses
R = 5.2  # distance between Jupiter and Sun in AU
omega = [0, 0, np.sqrt(G * (mass_sun + mass_jupiter) / R ** 3)]
orbital_period = 2 * np.pi / np.linalg.norm(omega)

# Take coordinate centre as centre of mass of system
r_s = np.array([-mass_jupiter * R / (mass_jupiter + mass_sun), 0, 0])  # Vector displacement from COM to Sun
r_j = np.array([mass_sun * R / (mass_jupiter + mass_sun), 0, 0])  # Vector displacement from COM to Jupiter


### TEST ONLY
# r_a = [2.5948051948051956, 4.50333209967908, 0]
# omega = [0, 0,0.5298767365821111]


def derivatives(t, y):
    r_a = y[:3]
    v_a = y[3:6]

    # Find displacement of asteroid from masses
    r_a_to_s = r_s - r_a
    r_a_to_j = r_j - r_a

    # Define force per unit asteroid mass due to gravity
    F = G * mass_sun * r_a_to_s / (np.linalg.norm(r_a_to_s) ** 3) + G * mass_jupiter * r_a_to_j / (
            np.linalg.norm(r_a_to_j) ** 3)

    # Define equations of motion
    r_a_dot = v_a
    v_a_dot = F - 2 * np.cross(omega, v_a) - np.cross(omega, np.cross(omega, r_a))

    return np.concatenate((r_a_dot, v_a_dot))


def plot_extras(ax1):
    # Plot Sun's location
    ax1.plot(r_s[0], r_s[1], 'yo')

    # Plot Jupiter's location
    ax1.plot(r_j[0], r_j[1], 'ro')

    # Plot COM
    ax1.plot(0, 0, 'b+')

    # Plot radius of jupiter
    orbit = plt.Circle((r_s[0], r_s[1]), R, fill=False, linewidth=1)
    ax1.add_artist(orbit)


def solve_orbit(initial_conditions, n_orbits):
    # initial_conditions is (r_a, v_a)
    t_max = 100 * orbital_period
    solution = scipy.integrate.solve_ivp(derivatives,
                                         t_span=(0, t_max),
                                         t_eval=np.linspace(0, t_max, 5000),
                                         y0=np.concatenate(initial_conditions))
    t, r_a, v_a = (solution.t, solution.y[:3], solution.y[3:6])
    return t, r_a, v_a


def plot_orbit():
    # Gives the default orbital view
    # Initial conditions of asteroid
    r_a_initial = [R * np.sin(np.pi / 6),
                   R * ((mass_sun - mass_jupiter) / (mass_sun + mass_jupiter)) * np.cos(np.pi / 6),
                   0]  # Asteroid vector displacement from COM
    r_a_initial = np.array(r_a_initial) + np.array([-0.05, +0.05, 0])
    print(r_a_initial)
    v_a_initial = [0, 0, 0]
    t, r_a, v_a = solve_orbit((r_a_initial, v_a_initial), 100)
    fig, ax1 = plt.subplots()
    ax1.set_xlim((-7, 7))
    ax1.set_ylim((-7, 7))
    ax1.plot(r_a[0], r_a[1])
    plot_extras(ax1)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()


def evaluate_wander():
    # Define a set of initial conditions in r_a, v_a phase space
    # Call solve orbit in parallel
    # Evaluate maximum wander using R_max
    grid_size = 2
    r_values = np.linspace(-0.05, +0.05, grid_size)
    r_lagrange_point = np.array([R * np.sin(np.pi / 6),
                                 R * ((mass_sun - mass_jupiter) / (mass_sun + mass_jupiter)) * np.cos(np.pi / 6),
                                 0])  # Asteroid vector displacement from COM
    # print((r_a_initial + r_lagrange_point))
    X, Y = np.meshgrid(r_values, r_values)
    r_max_array = np.zeros((grid_size, grid_size))
    for i in range(grid_size):
        for j in range(grid_size):
            r_a_initial = [X[i][j], Y[i][j], 0]
            v_a_initial = [0, 0, 0]
            t, r_a, v_a = solve_orbit((r_a_initial + r_lagrange_point, v_a_initial), 100)
            r_max = np.amax(np.sqrt((r_a[0] - r_lagrange_point[0]) ** 2 + (r_a[1] - r_lagrange_point[1]) ** 2))
            r_max_array[i][j] = r_max
    print(r_max_array)
    print(X)
    print(Y)
    # fig, ax1 = plt.subplots()
    # cp = ax1.map(X, Y, r_max_array)
    plt.show()


plot_orbit()
evaluate_wander()
