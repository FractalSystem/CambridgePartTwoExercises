import numpy.random as random
import numpy as np
import matplotlib
matplotlib.use('pdf') # required to compile on MCS over SSH
import matplotlib.pyplot as plt
import time

s = np.pi / 8


def integrand(input):
    return np.sin(np.sum(input)) * 10 ** 6


def estimate_average_f(N):
    # Uses Monte-Carlo routine to randomly sample f
    # Initiate pythonic list to store sampled values of the integrand, more efficient than np array here.
    f_list = []
    for i in range(N):
        # initiate 1D array of random values in range [0,s]
        inputs = random.rand(8) * s
        f = integrand(inputs)
        f_list.append(f)

    # Convert f_list to numpy array f_array for vectorisation
    f_array = np.array(f_list)
    f_squared_array = f_array ** 2
    average_f = (1 / N) * np.sum(f_array)
    average_f_squared = (1 / N) * np.sum(f_squared_array)
    return average_f, average_f_squared


def plot_errors(results_list):
    # Plot a graph of sampled error against N
    Ns = [dic.get("N") for dic in results_list]
    errors = [dic.get("sampled_error") for dic in results_list]
    fig, ax1 = plt.subplots()
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_ylabel("log(Sampled Error /No Units)")
    ax1.set_xlabel("log(N /No Units)")
    ax1.set_title("Plot of sampled error against N")
    ax1.plot(Ns, errors, label="Sampled error")

    # Calculate normalisation value to match expected relationship to results
    normalisation_k = errors[0] / Ns[0] ** -0.5
    ax1.plot(Ns, normalisation_k * np.power(Ns, -1 / 2), label="Theoretical relationship", linestyle="dashed")
    ax1.legend()
    plt.savefig("core_one_error_plot.pdf")


def plot_theoretical_error(results_list):
    Ns = [dic.get("N") for dic in results_list]
    errors = [dic.get("sampled_error") for dic in results_list]
    theoretical_errors = [dic.get("theoretical_error") for dic in results_list]
    fig, ax1 = plt.subplots()
    ax1.set_ylabel("Error /No Units")
    ax1.set_xlabel("N /No Units")
    ax1.plot(Ns, theoretical_errors, label="Theoretical error")
    ax1.plot(Ns, errors, label="Sampled error")
    ax1.set_title("Plot of theoretical and sampled error against N")
    ax1.legend()
    plt.savefig("supplementary_one_plot.pdf")


def main_loop(N_list, N_t=25):
    # Estimates the integral using Monte-Carlo techniques, by taking N samples, and repeating this process n_t times. This process is applied to each value N in N_list
    results_list = []
    volume = s ** 8

    for N in N_list:
        estimates_list = []
        theoretical_error_list = []
        for i in range(N_t):
            average_f, average_f_squared = estimate_average_f(N)
            estimate = average_f * volume
            theoretical_error = volume * ((average_f_squared - average_f ** 2) / N) ** (1 / 2)
            theoretical_error_list.append(theoretical_error)
            estimates_list.append(estimate)
        estimates_array = np.array(estimates_list)
        mean = np.mean(estimates_array)
        theoretical_error_mean = np.mean(theoretical_error_list)

        # RMS error in mean, from estimates
        sampled_error = np.std(estimates_array)
        results_list.append(
            {"N": N, "mean": mean, "sampled_error": sampled_error, "theoretical_error": theoretical_error_mean})
    return results_list


if __name__ == "__main__":
    N_list = np.arange(10, 400, 4, dtype=int).tolist()

    # Prove N^(-1/2) relationship
    print("Beginning error calculation.")
    results_list = main_loop(N_list, N_t=25)
    plot_errors(results_list)

    # Supplementary task one
    plot_theoretical_error(results_list)
    print("Finished error calculation.")

    # Estimate integral for N=20000, n_t=25
    N = 20000
    print(f"Estimating integral for N={N}.")
    t_start = time.time()
    results_list = main_loop([N], N_t=25)
    estimated_value = results_list[0].get("mean")
    sampled_error = results_list[0].get("sampled_error")
    print(
        f"Integral estimated to be {estimated_value} +/- {sampled_error} for N={N} in {round(time.time()-t_start, 3)}s.")
