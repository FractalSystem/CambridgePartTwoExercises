Core task 1:
The second order ODE was written as two coupled first order ODES and used scipy.integrate.solve_ivp to solve them
for q=F=0 and g=l. The code was tested by setting the initial theta to 0.01, allowing us to confirm the calculated
result agreed with the analytical result: theta = 0.01*cos(t). Energy was calculated for each time interval and
compared to the initial value over 10000 oscillations. The energy was found to not be conserved by the integrator,
and decreased exponentially. The results are seen in the file "rk4_energy_conservation_plot.png". The period was
also found for the range 0 < theta_initial < pi. For each initial theta value the period was predicted by averaging
the differences between the zero crossing points. The period was observed to increase from the 2*pi small angle
prediction rapidly and non-linearly. These results can be seen in the file "period_plot.png".

Core task 2:
The equation of motion was solved for q = 1,5,10 (keeping F=0)and the results plotted in the file
"varying_q_theta_plot.png". The behaviour corresponds to what we would expect, with q = 1 producing a light damping
profile and q=10 producing a heavily damped profile. The behaviour for F = 0.5, 1.2, 1.44, 1.465 under light damping
with q = 0.5 was plotted in the files "varying_F_theta(_dot)_plot.png". For high values of F the pendulum was
oberved to go "over-the-top" corresponding to values of |theta| > pi. The period was calculated numerically
from d(theta)/dt and found to increase as F increases, tending towards the driving period (=3*pi). The standard
deviation was also observed to increase as F increased, corresponding to the more chaotic motion observed.

Supplementary Task 1:
The equation of motion was solved for initial_theta = 0.2 and 0.20001. The solutions were found initially to agree very
closely but visibly diverged at around t=60s. The solutions then diverge rapidly as the initial_theta = 0.20001 goes
"over the top" while the other does not. This then affects the behaviour for all subsequent times. The code is appended
to the file "core_two.py" and results plotted in file "supplementary_one_varying_initial_theta.pdf".

Supplementary Task 2:
