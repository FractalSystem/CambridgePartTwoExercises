# Core functionality
######################################################

import math

import numpy as np
import scipy.integrate

import matplotlib.pyplot as plt

from time import sleep

######################################################
# user defined quantities

G = 4 * math.pi ** 2  # gravitational constant
ms = 1  # solar mass
mp = 0.001
R = 5.2  # astronomical units
T = R ** 1.5  # jupiter -> 11.8618 years

precision = 100  # evaluation points per orbit

######################################################
# derived quantities

rs = R * mp / (ms + mp)  # distance from origin to sun
rp = R * ms / (ms + mp)
w = 2 * math.pi / T  # angular velocity of frame
print(w)

# co-ordinates of lagrange point
lx = rp - R / 2
ly = math.sqrt(3) * R / 2


######################################################

# equations of motion
def derivs(t,y):
    rx, ry,rz, vx, vy, vz = y

    # print(y)
    # for convenience
    dp3 = ((rp - rx) ** 2 + ry ** 2) ** 1.5
    ds3 = ((rs + rx) ** 2 + ry ** 2) ** 1.5
    # print([-G * (ms * (rx + rs) / ds3 + mp * (rx - rp) / dp3), -G * (ms * ry / ds3 + mp * ry / dp3)])
    # print(rp)
    # print(rs)
    # print(rx)
    # print(ry)
    # print(dp3, ds3)
    # print(-G * (ms * (rx + rs) / ds3))
    # print(mp * (rx - rp) / dp3)
    # print(2 * np.linalg.norm(w) * vy)
    #
    # print(-G * (ms * (rx + rs) / ds3 + mp * (rx - rp) / dp3)+ 2 * w * vy )
    # print(-G * (ms * (rx + rs) / ds3 + mp * (rx - rp) / dp3) + 2 * w * vy + rx * w ** 2)
    # print(rx * w ** 2)
    # print(w)

    ret = (vx,
            vy, vz,
            -G * (ms * (rx + rs) / ds3 + mp * (rx - rp) / dp3) + 2 * w * vy + rx * w ** 2,
            -G * (ms * ry / ds3 + mp * ry / dp3) - 2 * w * vx + ry * w ** 2, 0)
    print(ret)
    # sleep(1000)
    return ret


# calculate a trajectory starting from y0
def orbit(y0, orbits):
    t_max = orbits * T
    print(t_max)
    print(np.array(y0))
    # y = scipy.integrate.odeint(derivs, y0, t)
    y  = scipy.integrate.solve_ivp(derivs,
                                     t_span=(0, t_max),
                                     t_eval=np.linspace(0, t_max, 5000),
                                     y0=np.array(y0))

    fig, ax1 = plt.subplots()
    ax1.set_xlim((-10, 10))
    ax1.set_ylim((-10, 10))

    ax1.plot(y.y[0], y.y[1])
    plt.show()

    return np.transpose(y)

orbit([2.5948051948051956, 4.50333209967908,0,0, 0,0], 100)
