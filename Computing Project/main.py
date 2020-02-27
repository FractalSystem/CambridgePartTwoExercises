import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

from time import sleep

G = 4 * np.pi ** 2
mass_jupiter = 0.001  # in solar masses
mass_sun = 1.0  # in solar masses
R = 5.2  # distance between Jupiter and Sun in AU
omega = [0, 0, np.sqrt(G * (mass_sun + mass_jupiter) / R ** 3)]
orbital_period = 2 * np.pi / np.linalg.norm(omega)

# take coordinate centre as centre of mass of system
r_s = np.array([-mass_jupiter * R / (mass_jupiter + mass_sun), 0, 0])  # Vector displacement from COM to Sun
r_j = np.array([mass_sun * R / (mass_jupiter + mass_sun), 0, 0])  # Vector displacement from COM to Jupiter

# Initial conditions of asteroid
r_a = [R * np.sin(np.pi / 6), R * ((mass_sun - mass_jupiter) / (mass_sun + mass_jupiter)) * np.cos(np.pi / 6),
       0]  # Asteroid vector displacement from COM

### TEST ONLY
# r_a = [2.5948051948051956, 4.50333209967908, 0]
# omega = [0, 0,0.5298767365821111]

v_a = [0, 0, 0]


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


t_max = 100 * orbital_period
print(t_max)
print(np.concatenate((r_a, v_a)))
solution = scipy.integrate.solve_ivp(derivatives,
                                     t_span=(0, t_max),
                                     t_eval=np.linspace(0, t_max, 5000),
                                     y0=np.concatenate((r_a, v_a)))
t, r_a, v_a = (solution.t, solution.y[:3], solution.y[3:6])
fig, ax1 = plt.subplots()
ax1.set_xlim((-10, 10))
ax1.set_ylim((-10, 10))
ax1.plot(r_a[0], r_a[1])
plt.show()
