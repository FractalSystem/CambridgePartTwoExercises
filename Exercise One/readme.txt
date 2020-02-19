Core task 1:
The integral was estimated using Monte-Carlo techniques which were implemented in a vectorised form. The error in the
integral was estimated from the standard deviation of 25 estimates for each N. This sampled error was confirmed to
be proportional to N^(-1/2). The results are plotted in "error_plot.pdf".
For N=20000, the integral was estimated to be 537.1 +/- 0.3 in 8.2s over 25 iterations, which was deemed a reasonable time.

Core task 2:
The fresnel integrals were evaluated using the scipy.integrate.quad() function for values of u where -4*pi < u < 4*pi.
A Cornu spiral was plotted and saved in "cornu_spiral_plot.pdf".

Supplementary task 1:
Theoretical error and sampled error were observed to converge quickly as N increased, and to be in agreement for large N.
Source code is incorporated into core_one.py. Results plotted in "supplementary_one_plot.pdf"

Supplementary task 2:
Source in file "supplementary_two.py". Results plotted in "supplementary_two_plot.pdf"