import numpy as np
import scipy.integrate
import matplotlib

matplotlib.use('agg')
import matplotlib.pyplot as plt


def plot_q(solution_dic_list):
    solution_list = [d.get("solution") for d in solution_dic_list]
    q_list = [d.get("q") for d in solution_dic_list]
    fig, (ax1, ax2) = plt.subplots(2, figsize=[6.4, 9.6])
    for i, q in enumerate(q_list):
        solution = solution_list[i]
        t, theta, theta_dot = (solution.t, solution.y[0], solution.y[1])
        ax1.plot(t, theta, label=f"q = {str(q)}")
        ax2.plot(t, theta_dot, label=f"q = {str(q)}")
    ax1.set_xlabel("Time /s")
    ax1.set_ylabel("Theta /rad")
    ax1.set_title("Plot of Theta Against Time for Three Values of q. F = 0")
    ax1.legend()
    ax2.set_xlabel("Time /s")
    ax2.set_ylabel("d(Theta)/dt /rad/s")
    ax2.set_title("Plot of d(Theta)/dt Against Time for Three Values of q. F = 0")
    ax2.legend()
    plt.tight_layout()  # prevents cut off of y label
    plt.savefig("core_two_varying_q_plot.pdf")


def plot_f_theta(solution_dic_list):
    solution_list = [d.get("solution") for d in solution_dic_list]
    f_list = [d.get("f") for d in solution_dic_list]
    fig = plt.figure(figsize=[6.4 * 4, 4.8 * 2 + 1])
    (ax_theta, ax_theta_dot) = fig.subplots(2, 4)
    for i, f in enumerate(f_list):
        solution = solution_list[i]
        t, theta, theta_dot = (solution.t, solution.y[0], solution.y[1])
        ax1 = ax_theta[i]
        ax2 = ax_theta_dot[i]
        ax1.plot(t, theta, label=f"F = {str(f)}")
        ax2.plot(t, theta_dot, label=f"F = {str(f)}")

        # Graph bits
        ax1.set_xlabel("Time /s")
        ax1.set_ylabel("Theta /rad")
        ax1.set_ylim([-40, 10])
        ax1.legend()

        ax2.set_xlabel("Time /s")
        ax2.set_ylabel("d(Theta)/dt /rad/s")
        ax2.set_ylim([-2.7, 2.7])
        ax2.legend()
    # plt.savefig("varying_F_theta_plot.png")
    fig.suptitle("Plots of theta and d(Theta)/dt against time for four values of F. q = 0.5")
    plt.tight_layout()  # prevents cut off of y label
    plt.subplots_adjust(top=0.95)
    plt.savefig("core_two_varying_F_plot.pdf")


def plot_supp_one(solution_dic_list):
    solution_list = [d.get("solution") for d in solution_dic_list]
    initial_theta_list = [d.get("initial_theta") for d in solution_dic_list]
    fig = plt.figure(figsize=[6.4, 9.6])
    (ax1, ax2) = fig.subplots(2)
    for i, initial_theta in enumerate(initial_theta_list):
        solution = solution_list[i]
        t, theta, theta_dot = (solution.t, solution.y[0], solution.y[1])
        ax1.plot(t, theta, label=f"Initial theta = {str(initial_theta)} rad")
        ax2.plot(t, theta_dot, label=f"Initial theta = {str(initial_theta)} rad")
    ax1.set_xlabel("Time /s")
    ax1.set_ylabel("Theta /rad")
    ax1.set_title("Plot of Theta Against Time for Varying Initial Theta")
    ax1.legend()
    ax2.set_xlabel("Time /s")
    ax2.set_ylabel("d(Theta)/dt /rad/s")
    ax2.set_title("Plot of d(Theta)/dt Against Time for Varying Initial Theta")
    ax2.legend()
    plt.tight_layout()  # prevents cut off of y label
    plt.savefig("supplementary_one_varying_initial_theta.pdf")


def derivatives(t, y, q, F):
    return [y[1], -np.sin(y[0]) - q * y[1] + F * np.sin(2 / 3 * t)]


def solve(y0, t_max, q, F):
    """
    function to solve the ODE with initial conditions and constants as the arguments
    """
    solution = scipy.integrate.solve_ivp(
        fun=derivatives,
        t_span=(0, t_max),
        y0=y0,
        args=(q, F,),
        t_eval=np.linspace(0, t_max, 10000),
    )
    # for i in range(len(solution.t)):
    #     print(f"{solution.t[i]}, {solution.y[0, i]}, {solution.y[1,i]}")
    return solution


def get_period(solution):
    crossing_times = []
    negative_area = False
    for i in range(len(solution.t)):
        if not negative_area:
            if solution.y[1, i] < 0:
                # Zero crossing
                crossing_time = solution.t[i]
                crossing_times.append(crossing_time)
                negative_area = True
        else:
            if solution.y[1, i] > 0:
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
    std_dev = np.std(periods)
    return mean_period, std_dev


def supplementary_two(q_values, f_values):
    assert len(q_values) > 1
    assert len(f_values) > 1
    # fig = plt.figure(figsize=[6.4*2, 9.6*2])
    fig, axs = plt.subplots(len(f_values), len(q_values), figsize=[3.2 * len(q_values), 2.8 * len(f_values)])
    for i, f in enumerate(f_values):
        for j, q in enumerate(q_values):
            # initial theta is 0.01
            solution = solve((0.01, 0), 10 * 2 * np.pi, q, f)
            t, theta, theta_dot = (solution.t, solution.y[0], solution.y[1])
            # ax[row, col]
            ax = axs[i, j]
            if j == 0:
                axs[i, 0].set_ylabel(f"F={f}\n\nAngular speed (rad/s)")
            if i == len(f_values) - 1:
                axs[i, j].set_xlabel(f"Angle (rad)\n\nq={q}")
            ax.plot(theta, theta_dot)
    fig.suptitle("Plot of angular speed against angle for various values of q and F\n\n(Note that x and y axis labels apply to all plots)")
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.2, top=0.92, bottom=0.1)
    plt.savefig("supplementary_two.pdf")
    # plt.show()


if __name__ == "__main__":
    print("Plotting results for varying q.")
    q_list = [1, 5, 10]
    solution_list = []
    for q in q_list:
        solution = solve((0.01, 0), 10 * 2 * np.pi, q, 0)
        solution_list.append({"solution": solution, "q": q})
    plot_q(solution_list)

    print("Plotting results for varying F.")
    f_list = [0.5, 1.2, 1.44, 1.465]
    period_list = []
    solution_dic_list = []
    for f in f_list:
        solution = solve((0.01, 0), 10 * 2 * np.pi, 0.5, f)
        mean_period, std_dev = get_period(solution)
        period_list.append((mean_period, std_dev))
        solution_dic_list.append({"solution": solution, "f": f})
    plot_f_theta(solution_dic_list)
    print(f"[(period, std_dev)] : {period_list}")

    # Supplementary exercise one
    print("Plotting results for varying initial_theta.")
    initial_theta_list = [0.2, 0.20001]
    solution_dic_list = []
    for initial_theta in initial_theta_list:
        solution = solve((initial_theta, 0), 30 * 2 * np.pi, 0.5, 1.2)
        solution_dic_list.append({"solution": solution, "initial_theta": initial_theta})
    plot_supp_one(solution_dic_list)

    # # Supplementary exercise two
    print("Plotting results for supplementary exercise two.")
    q_values = [0, 0.5, 5.0]
    f_values = [0, 0.5, 1.0, 1.44, 2]
    # q_values = [0,1]
    # f_values = [0,1,2]
    supplementary_two(q_values, f_values)
