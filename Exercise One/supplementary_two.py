import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

wavelength = 1*10**-2
slit_width = 10*10**-2


def function_c(x, distance):
    return np.cos((np.pi*x**2)/(wavelength*distance))

def function_s(x, distance):
    return np.sin((np.pi*x**2)/(wavelength*distance))


def evaluate_pattern(distance):
    x1_array = np.linspace(0, -slit_width, 1000) # distance from right edge
    x2_array = np.linspace(slit_width, 0, 1000) # distance from left edge
    c_list = []
    s_list = []
    for x1,x2 in zip(x1_array, x2_array):
        c, err = integrate.quad(function_c, x1 ,x2, args=(distance,))
        s, err = integrate.quad(function_s, x1 ,x2, args=(distance,))
        c_list.append(c)
        s_list.append(s)
    c_array = np.array(c_list)
    s_array = np.array(s_list)
    amplitude_array = (c_array**2+s_array**2)
    phase_array = np.tan(s_array/c_array)
    return amplitude_array, phase_array
    # plt.plot(x1_array, amplitude)
    # plt.show()
    # plot(c_list, s_list)

if __name__ == "__main__":
    distances = [30*10**-2, 50*10**-2, 100*10**-2]
    x1_array = np.linspace(0, -slit_width, 1000)
    fig = plt.figure(figsize= [6.4, 9.6])
    (ax1,ax2) = fig.subplots(2)
    for distance in distances:
        amplitude_array, phase_array = evaluate_pattern(distance)
        ax1.plot(-x1_array ,amplitude_array, label=f"D = {int(distance*10**2)}cm")
        ax2.plot(-x1_array, phase_array, label=f"D = {int(distance*10**2)}cm")

    ax1.legend()
    ax1.set_xlabel("Distance from edge of slit /cm")
    ax1.set_ylabel("Relative amplitude")
    ax1.set_title("Plot of relative amplitude across the screen")
    ax2.legend()
    ax2.set_xlabel("Distance from edge of slit /cm")
    ax2.set_ylabel("Relative phase")
    ax2.set_title("Plot of the relative phase across the screen")

    plt.savefig("supplementary_two.pdf")

