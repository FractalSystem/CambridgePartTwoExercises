import numpy as np
import matplotlib
matplotlib.use('pdf') # required to compile on MCS over SSH
import matplotlib.pyplot as plt
from scipy import integrate


def function_c(x):
    return np.cos((np.pi*x**2)/2)

def function_s(x):
    return np.sin((np.pi*x**2)/2)

def plot(c_list, s_list):
    fig, ax = plt.subplots()
    ax.set_xlabel("C(u)")
    ax.set_ylabel("S(u)")
    ax.set_title("Plot of the Cornu spiral")
    plt.plot(c_list, s_list)
    plt.plot(0.5, 0.5, "r+")
    plt.plot(-0.5, -0.5, "r+")
    plt.savefig("core_two_cornu_spiral_plot.pdf")


def main():
    # u_array is an array of u values to calculate c and s for
    u_array = np.linspace(-4*np.pi, 4*np.pi, 2000)
    c_list = []
    s_list = []
    for u in u_array:
        c, err = integrate.quad(function_c, 0 ,u)
        s, err = integrate.quad(function_s, 0 ,u)
        c_list.append(c)
        s_list.append(s)
    plot(c_list, s_list)

if __name__ == "__main__":
    main()