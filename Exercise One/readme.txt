Core task 1 ("core_one.py"):
The integral was estimated using Monte-Carlo techniques which were implemented in a vectorised form. The error in the
integral was estimated from the standard deviation of 25 estimates for each N. This sampled error was confirmed to
be proportional to N^(-1/2). The results are plotted in "core_one_error_plot.pdf".
For N=20000, the integral was estimated to be 537.1 +/- 0.3 in 8.2s over 25 iterations, which was deemed a reasonable time.

Core task 2 ("core_two.py"):
The fresnel integrals were evaluated using the scipy.integrate.quad() function for values of u where -4*pi < u < 4*pi.
A Cornu spiral was plotted and saved in "core_two_cornu_spiral_plot.pdf". The Cornu spiral converges on the points
(0.5, 0.5) and (-0.5, -0.5) as expected.

Supplementary task 1 ("core_one.py"):
The theoretical error was calculated for each N using the formula given in the handout. Theoretical error and sampled
error were observed to agree closely. Sampled error was seen to fluctuate more than theoretical error, as it is more
susceptible to variations in the randomly selected initial values than the theoretical error.
Source code is incorporated into "core_one.py". Results plotted in "supplementary_one_plot.pdf"

Supplementary task 2 ("supplementary_two.py"):
The near field complex amplitude integral was evaluated for three different values of distance, D. Relative amplitude
and phase were plotted as functions of distance from centre of screen. The results correspond with what we would expect,
with the peaks in amplitude growing further apart as D increases.
Source in file "supplementary_two.py". Results plotted in "supplementary_two_plot.pdf"