import numpy as np
import matplotlib
matplotlib.use('pdf') # required to compile on MCS over SSH
import matplotlib.pyplot as plt
from scipy import integrate

wavelength = 1
slit_width = 10


def function_c(x, distance):
    return np.cos((np.pi * (x ** 2)) / (wavelength * distance))


def function_s(x, distance):
    return np.sin((np.pi * (x ** 2)) / (wavelength * distance))


def evaluate_pattern(distance):
    x1_array = np.linspace(-30 - slit_width / 2, +30 - slit_width / 2, 1000)  # distance from right edge
    x2_array = x1_array + slit_width  # distance from left edge
    c_list = []
    s_list = []
    for x1, x2 in zip(x1_array, x2_array):
        c, err = integrate.quad(function_c, x1, x2, args=(distance,))
        s, err = integrate.quad(function_s, x1, x2, args=(distance,))
        c_list.append(c)
        s_list.append(s)
    c_array = np.array(c_list)
    s_array = np.array(s_list)
    amplitude_array = (c_array ** 2 + s_array ** 2) ** 0.5
    amplitude_array = amplitude_array / np.max(amplitude_array)
    phase_array = np.arctan(s_array / c_array)
    return amplitude_array, phase_array


if __name__ == "__main__":
    distances = [30, 50, 100]
    x1_array = np.linspace(-30, +30, 1000)
    fig = plt.figure(figsize=[3.2 * 6, 2.8 * 4])
    (ax_ampl, ax_phase) = fig.subplots(2, 3)
    for i, distance in enumerate(distances):
        amplitude_array, phase_array = evaluate_pattern(distance)
        ax_ampl[i].plot(-x1_array, amplitude_array, label=f"D = {int(distance)}cm")
        ax_ampl[i].set_xlabel("Distance from screen centre /cm")
        ax_ampl[i].set_ylabel("Relative amplitude /No Units")
        ax_ampl[i].legend()
        ax_phase[i].plot(-x1_array, phase_array, label=f"D = {int(distance)}cm")
        ax_phase[i].set_xlabel("Distance from screen centre /cm")
        ax_phase[i].set_ylabel("Relative phase /rad")
        ax_phase[i].legend()

    fig.suptitle("Plots of relative amplitude and phase against distance from centre of screen for three values of D")
    plt.savefig("supplementary_two_plot.pdf")
