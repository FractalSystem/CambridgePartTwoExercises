import numpy as np
import scipy.integrate
import matplotlib

matplotlib.use('agg')  # required to compile on MCS over SSH
import matplotlib.pyplot as plt


def small_angle_plot(solution):
    """
    Function to plot small angle approx
    """
    t, theta, theta_dot = (solution.t, solution.y[0], solution.y[1])
    fig, ax1 = plt.subplots()
    ax1.plot(t, theta, label="Theta (rad)")
    ax1.plot(t, 0.01 * np.cos(t), label="Theoretical result", linestyle="dashed")
    ax1.set_xlabel("Time /s")
    ax1.set_ylabel("Theta /rad")
    ax1.set_title("Plot of theta against time for q = F = 0")
    # ax1.set_xlim(left = 998*2*np.pi, right = 1000*2*np.pi)
    ax1.legend()
    plt.tight_layout()  # prevents cut off of y label
    plt.savefig("core_one_small_angle_plot.pdf")


def plot_energy(solution):
    """
    Function to plot energy evolution over time. This tests how well the integrator coserves energy
    """
    t, theta, theta_dot = (solution.t, solution.y[0], solution.y[1])
    # energy is sum of KE+PE, l=g=10
    g = 10
    l = g
    m = 1
    energy = 1 / 2 * m * (theta_dot * l) ** 2 + m * g * (l - l * np.cos(theta))
    initial_energy = 1 / 2 * m * (0 * l) ** 2 + m * g * (l - l * np.cos(0.01))

    fig, (ax1, ax2) = plt.subplots(2, figsize=[3.2 * 2, 2.8 * 4])

    # Plot energy calculated for each element of the solution
    ax1.plot(t, energy, label="RK4 Energy")

    # Energy is conserved
    ax1.plot(t, [initial_energy] * (len(t)), label="Theoretical energy")
    ax1.legend()
    ax1.set_xlabel("Time /s")
    ax1.set_ylabel("Energy /J/kg")
    ax1.set_title("Plot of energy against time for theoretical and \nRK4 solutions for 10,000 oscillations")
    ax2.plot(t, (1 - energy / initial_energy) * 100)
    ax2.set_xlabel("Time /s")
    ax2.set_ylabel("Percentage deviation /%")
    ax2.set_title("Plot of percentage deviation of RK4 energy from \ntheoretical energy for 10,000 oscillations")
    plt.tight_layout()  # prevents cut off of y label
    plt.savefig("core_one_energy_conservation.pdf")


def plot_period_theta(theta_initials, periods):
    fig, ax1 = plt.subplots()
    ax1.plot(theta_initials, periods, label="Period")
    ax1.set_xlabel("Initial theta /rad")
    ax1.set_ylabel("Period /s")
    ax1.set_title("Plot of Initial Theta Against Period")
    plt.tight_layout()  # prevents cut off of y label
    plt.savefig("core_one_period_plot.pdf")


def derivatives(t, y, q, F):
    return [y[1], -np.sin(y[0]) - q * y[1] + F * np.sin(2 / 3 * t)]


def solve(y0, t_max, q, F):
    """
    function to solve the ODE with initial conditions and constants as the arguments. t_max defines the time range to solve over.
    """
    solution = scipy.integrate.solve_ivp(
        fun=derivatives,
        t_span=(0, t_max),
        y0=y0,
        args=(q, F,),
        t_eval=np.linspace(0, t_max, 10000),
    )
    return solution


def solve_period(theta_initial):
    """
    Find period for given initial theta
    """
    y0 = (theta_initial, 0)
    t_max = 100 * 2 * np.pi
    q = 0
    F = 0
    solution = scipy.integrate.solve_ivp(
        fun=derivatives,
        t_span=(0, t_max),
        y0=y0,
        args=(q, F,),
        t_eval=np.linspace(0, t_max, 10000),
        rtol=1e-6,
        atol=1e-6
    )

    # routine to find mean period by sampling all zero crossing points of theta and averaging the differences.
    crossing_times = []
    negative_area = False
    for i in range(len(solution.t)):
        if not negative_area:
            if solution.y[0, i] < 0:
                # Zero crossing
                crossing_time = solution.t[i]
                crossing_times.append(crossing_time)
                negative_area = True
        else:
            if solution.y[0, i] > 0:
                # Zero crossing
                crossing_time = solution.t[i]
                crossing_times.append(crossing_time)
                negative_area = False
    periods = []
    previous_crossing = crossing_times[0]
    for i in range(len(crossing_times) - 1):
        periods.append((crossing_times[i + 1] - previous_crossing) * 2)
        previous_crossing = crossing_times[i + 1]
    mean_period = np.mean(periods)
    return mean_period


def period_loop():
    """
    Call solve_period for a range of initial thetas and return periods
    """
    theta_initials = np.linspace(0.01, np.pi - 0.01, 64)
    periods = []
    for theta_initial in theta_initials:
        period = solve_period(theta_initial)
        periods.append(period)
    # for i in range(len(theta_initials)):
    #     print(f"{theta_initials[i]}, {periods[i]}")
    return theta_initials, periods


if __name__ == "__main__":
    # Small angle approximation
    print("Testing small angle solution.")
    solution = solve((0.01, 0), 4 * 2 * np.pi, 0, 0)
    small_angle_plot(solution)

    # Conservation of energy
    print("Testing conservation of energy.")
    solution = solve((0.01, 0), 10000 * 2 * np.pi, 0, 0)
    plot_energy(solution)

    # Evaluate period for a range of initial thetas
    print("Plotting period against initial theta.")
    theta_initials, periods = period_loop()
    plot_period_theta(theta_initials, periods)

    # find period for theta initial = pi/2 as requested
    print(f"Period for initial theta = pi/2: {solve_period(np.pi / 2)}s")  # answer = 7.411925566125717 s
