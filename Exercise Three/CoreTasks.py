import numpy as np
import matplotlib.pyplot as plt

# Define constants
N = 2 ** 17  # Number of sample points in aperture. 2^n most efficient for FFT.


def plot_fft_core_one(near_field, Y, y, theoretical_Y, theoretical_y):
    """
    Plotting routine for Core Task 1.
    :param near_field: (bool) Flag stating whether we are plotting in the near field regime.
    :param Y: (numpy.ndarray) Naturally ordered FFT of the aperture function.
    :param y: (numpy.ndarray) Naturally ordered y-space coordinates, stretched to be have correct units (mm).
    :param theoretical_Y: (numpy.ndarray) Theoretical wavefunction for single slit array (sinc function).
    :param theoretical_y: (numpy.ndarray) y-coordinates of samples of theoretical wavefunction.
    :return:
    """

    fig, ax1 = plt.subplots()

    # Normalise wavefunctions and square magnitude to get intensity
    Y_intensity = (abs(Y) / max(abs(Y))) ** 2
    theoretical_Y_intensity = (theoretical_Y / max(theoretical_Y)) ** 2

    # Plot over limited central range (to remove boring low intensity edges). Set partial transparency of topmost line.
    plot_range = 300  # Number of points to plot either side of centre
    ax1.plot(y[int(N / 2) - plot_range:int(N / 2) + plot_range],
             Y_intensity[int(N / 2) - plot_range:int(N / 2) + plot_range], label="FFT")
    if not near_field:
        ax1.plot(theoretical_y[int(N / 2) - plot_range:int(N / 2) + plot_range],
                 theoretical_Y_intensity[int(N / 2) - plot_range:int(N / 2) + plot_range], label="Theoretical",
                 alpha=0.7)

    # Set plot properties
    ax1.set_xlabel("y /mm")
    ax1.set_ylabel("Normalised intensity /No units")
    if not near_field:
        ax1.set_title("Plot of intensity against y for a single slit in the far-field regime")
        ax1.legend()
    else:
        ax1.set_title("Plot of intensity against y for a single slit in the near-field regime")
        ax1.legend()

    # Save plot to pdf
    if near_field:
        filename = "fresnel_supp_core_one_plot.pdf"
    else:
        filename = "core_one_plot.pdf"
    plt.tight_layout()  # prevents cut off of y label
    plt.savefig(filename)


def plot_FFT_core_two(near_field, y, Y):
    """
    Plotting routine for Core Task 2.
    :param near_field: (bool) Flag stating whether we are plotting in the near field regime.
    :param Y: (numpy.ndarray) Naturally ordered FFT of the aperture function.
    :param y: (numpy.ndarray) Naturally ordered y-space coordinates, stretched to be have correct units (mm).
    :return:
    """
    fig, ax1 = plt.subplots()
    Y_intensity = (abs(Y) / max(abs(Y))) ** 2  # Normalise wavefunction and square magnitude

    # Plot over limited central range (to remove boring low intensity edges)
    plot_range = 800  # Number of points to plot either side of centre
    ax1.plot(y[int(N / 2) - plot_range:int(N / 2) + plot_range],
             Y_intensity[int(N / 2) - plot_range:int(N / 2) + plot_range])

    # Set plot properties
    ax1.set_xlabel("y /mm")
    ax1.set_ylabel("Normalised intensity /No units")
    if not near_field:
        ax1.set_title("Plot of intensity against y for a grating in the far-field regime")
    else:
        ax1.set_title("Plot of intensity against y for a grating in the near-field regime")

    # Save plot to pdf
    if near_field:
        filename = "fresnel_supp_core_two_plot.pdf"
    else:
        filename = "core_two_plot.pdf"
    plt.tight_layout()  # prevents cut off of y label
    plt.savefig(filename)


def generate_slit_aperture(L, d):
    """
    Function generates a single slit of width d in the middle of an aperture of extent L.
    The slit will be at the centre of the aperture array, and will have transmission amplitude = 1.0.

    :param L: (float) Extent (mm)
    :param d: (float) Slit width (mm)
    :return: (numpy.ndarray) Aperture array
    """

    delta = L / N

    # Initialise zero valued aperture
    aperture = np.zeros(N, dtype=np.complex128)

    # Set slit between -d/2, +d/2
    d_min_index = int(round((-d / 2 + L / 2) / delta))
    d_max_index = int(round((d / 2 + L / 2) / delta))
    aperture[d_min_index:d_max_index] = 1.0
    return aperture


def generate_sin_aperture(L, d, m, s):
    """
    Function generates sinusoidal aperture of phase (m/2)sin(2*pi*x/s) and transmission amplitude 1.0. The overall
    aperture function is then of the form 1.0*exp(i*phase) where i is sqrt(-1). The grating will be generated in the
    centre of the aperture.
    :param L: (float) Extent (mm)
    :param d: (float) Slit width (mm)
    :param m: (float) Parameter in grating phase (dimensionless)
    :param s: (float) Parameter in grating phase (dimensionless)
    :return:
    """
    delta = L / N

    # Calculate grating function across entire aperture
    x = np.arange(-L / 2, L / 2, delta)
    grating = 1.0 * np.exp((m / 2) * np.sin(2 * np.pi * x / s) * 1j)
    # grating = (m / 2) * np.sin(2 * np.pi * x / s)

    # Initialise zero valued aperture
    aperture = np.zeros(N, dtype=np.complex128)

    # Set aperture equal to grating function between -d/2, +d/2
    d_min_index = int(round((-d / 2 + L / 2) / delta))
    d_max_index = int(round((d / 2 + L / 2) / delta))
    aperture[d_min_index:d_max_index] = grating[d_min_index:d_max_index]
    return aperture


def core_one(near_field=False):
    """
    Main program routine for Core Task 1. Defines problem parameters and calculates far-field diffraction pattern for
    single slit aperture.
    :param near_field: (bool) Flag stating whether we are plotting in the near field regime.
    :return:
    """

    # Define problem parameters
    L = 5.0  # Aperture extent in mm
    d = 0.1  # Slit width in mm
    D = 1000  # Screen distance in mm
    wavelength = 500 * 10 ** (-6)  # Wavelength in mm
    delta = L / N

    # Generate aperture
    aperture = generate_slit_aperture(L, d)

    # Modify aperture for near field if near_field flag is true
    if near_field:
        # Redefine parameter
        D = 5  # Screen distance in mm

        # Modify aperture function for near field
        k = 2 * np.pi / wavelength
        x = np.arange(-L / 2, L / 2, delta)
        aperture = aperture * np.exp(1j * k * x ** 2 / (2 * D))

    # Perform FFT on aperture
    Y = np.fft.fft(aperture)  # far field A=A'
    f = np.fft.fftfreq(N, delta)  # gives spaced values of j corresponding to coords

    # Convert f into y (proper coordinates)
    y = (-(D * wavelength)) * f

    # Swap spectrum halfs to create naturally ordered arrays (avoids x=0 line).
    Y = np.fft.fftshift(Y)
    y = np.fft.fftshift(y)

    # Theoretical result for a single slit
    theoretical_y = np.linspace(min(y), max(y), N)  # Sample theoretical solution only in range of numerical solution.
    q = (2 * np.pi / wavelength) * np.sin(np.arctan(theoretical_y / D))
    theoretical_Y = d * np.sin(q * d / 2) / (q * d / 2)

    # Plot numerical and theoretical solutions
    plot_fft_core_one(near_field, Y, y, theoretical_Y, theoretical_y)


def core_two(near_field=False):
    """
    Main routine for Core Task 2. Defines problem parameters and calculates far-field diffraction pattern for
    sinusoidal aperture.
    :param near_field: (bool) Flag stating whether we are plotting in the near field regime.
    :return:
    """

    # Define parameters
    L = 10.0  # aperture extent in mm
    d = 2.0  # slit width in mm
    m = 8.0  # Aperture phase prefactor in mm
    s = 0.1  # spacing of phase maxima in mm
    fresnel_distance = 8000.0  # d^2/wavelength in mm
    D = 10000  # screen distance in mm
    delta = L / N
    wavelength = d ** 2 / fresnel_distance  # wavelength in mm

    # Generate aperture
    aperture = generate_sin_aperture(L, d, m, s)

    # Modify aperture for near field if near_field flag is true
    if near_field:
        # Redefine parameter
        D = 500  # screen distance in mm

        # Modify aperture function
        k = 2 * np.pi / wavelength
        x = np.arange(-L / 2, L / 2, delta)
        aperture = aperture * np.exp(1j * k * x ** 2 / (2 * D))

    # Perform FFT on aperture
    Y = np.fft.fft(aperture)  # far field A=A'
    f = np.fft.fftfreq(N, delta)  # gives spaced values of j corresponding to coords

    # Convert f into y (proper coordinates)
    y = (-(D * wavelength)) * f

    # Swap spectrum halfs to create naturally ordered arrays (avoids x=0 line).
    Y = np.fft.fftshift(Y)
    y = np.fft.fftshift(y)

    # Plot solution
    plot_FFT_core_two(near_field, y, Y)


if __name__ == "__main__":
    # Run core tasks
    core_one()
    core_two()

    # Run supplementary task for Fresnel regime
    core_one(near_field=True)
    core_two(near_field=True)
