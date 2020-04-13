import numpy as np

# Define constants
G = 4 * np.pi ** 2
MASS_JUPITER = 0.001  # in solar masses
MASS_SUN = 1.0  # in solar masses
R = 5.2  # distance between Jupiter and Sun in AU

# ffmpeg is required to save animations. Under Windows the exe can be downloaded from https://ffmpeg.zeranoe.com/builds/
FFMPEG_PATH = r"C:\ffmpeg.exe"

# print(MASS_JUPITER*MASS_SUN/(MASS_JUPITER+MASS_SUN))
# omega = np.sqrt(G * (MASS_SUN + MASS_JUPITER) / R ** 3)
# print(2*np.pi/omega)