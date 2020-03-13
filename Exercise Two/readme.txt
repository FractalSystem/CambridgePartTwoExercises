Core task 1 ("core_one.py"):
The second order ODE was written as two coupled first order ODES and scipy.integrate.solve_ivp was used to solve them
for q = F = 0 and g = l. The code was tested by setting the initial theta to 0.01, allowing us to confirm the calculated
result agreed with the small-angle analytical result: theta = 0.01*cos(t) for small times. A plot of the same
calculation at large t shows that the solutions do not agree, although adjusting the atol and rtol parameters of
solve_ivp can improve the situation at the expense of longer computation time. These graphs are plotted in the files
"core_one_small_angle_plot(_large_t).pdf".
Energy was calculated for each time interval and compared to the initial value over 10000 oscillations. The energy was
found to not be conserved by the integrator, but instead decreased. The results are plotted in the file
"rk4_energy_conservation_plot.pdf".
The period was also found for the range 0 < theta_initial < pi. For each initial theta value the period was predicted by
averaging the differences between the zero crossing points. The period was observed to increase from the 2*pi small
angle prediction non-linearly. These results can be seen in the file "core_one_period_plot.pdf".
The period of theta_initial = pi/2 was calculated to be 7.41s.

Core task 2 ("core_two.py"):
The equation of motion was solved for q = 1,5,10 (keeping F=0) and the results plotted in the file
"core_two_varying_q_plot.pdf". The behaviour corresponds to what we would expect, with q = 1 producing a lightly
damped profile and q=10 producing a heavily damped profile.
The behaviour for F = 0.5, 1.2, 1.44, 1.465 under light damping with q = 0.5 was plotted in the file
"core_two_varying_F_plot.pdf". For F = 1.2, 1.44, 1.465 values of F the pendulum was observed to go "over-the-top"
corresponding to values of |theta| > pi. The period was calculated numerically from d(theta)/dt, although it has little
meaning for F = 1.2, 1.44, 1.465 as the motion is not periodic. For F = 0.5 the pendulum stabilises with a much larger
theta than in the undamped case, as expected.

Supplementary Task 1 ("core_two.py"):
The equation of motion was solved for initial_theta = 0.2 and 0.20001. The solutions were found initially to agree very
closely but visibly diverged at around t=60s. The solutions then diverge rapidly as the initial_theta = 0.20001 goes
"over the top" while the other does not. This then affects the behaviour for all subsequent times. The code is appended
to the file "core_two.py" and results are plotted in file "supplementary_one_varying_initial_theta.pdf".

Supplementary Task 2 ("core_two.py"):
Plots of angular speed against angle (theta) were plotted for a range of q and F values. This revealed a number of
interesting regimes. The plot can be found in the file "supplementary_two.pdf".
The F=q=0 graph is circular as this represents the pi/2 phase difference between angle and angular speed.
The plots become less chaotic as q increases, as damping suppresses the effect of the driving force, requiring a higher
F to produce a chaotic effect. In the case of q = 5 the pendulum is insensitive to all values of F covered by the plots,
producing approximately periodic motion as a result.
For q = 0.5 the motion becomes chaotic for increasing F, starting with F = 1.44. Note that the graph for F = 1.44,
q = 0.5 goes backwards (to negative theta) while several others go forwards (to positive theta). This shows how
sensitive the pendulum is to the initial conditions.