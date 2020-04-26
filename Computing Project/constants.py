#######################
# File defining all constant values. Imported by every file in the program to ensure consistency of constant values.
#######################

import numpy as np

# Define constants
G = 4 * np.pi ** 2
MASS_JUPITER = 0.001  # in solar masses
MASS_SUN = 1.0  # in solar masses
R = 5.2  # distance between Jupiter and Sun in AU

# ffmpeg is required to save animations. Under Windows the exe can be downloaded from https://ffmpeg.zeranoe.com/builds/
FFMPEG_PATH = r"C:\ffmpeg.exe"
