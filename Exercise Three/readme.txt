Core Task 1 ("CoreTasks.py"):
The diffraction pattern for a one dimensional slit of width d in the centre of an aperture of extent L was calculated
by taking the FFT (fast fourier transform) of the sampled aperture function. This calculation was performed assuming the
Fraunhofer limit, which is reasonable as d >> x_max**2/wavelength. The theoretical intensity was also plotted
and found to agree very closely with the calculated value.

Core Task 2 ("CoreTasks.py"):
The aperture function of a sinusoidal phase grating was calculated, resulting in a complex function. As above, the
Fraunhofer diffraction pattern was calculated using a FFT technique (although the Fraunhofer limit was not strictly
satisfied, as 10 >> 8 is not true). The resulting pattern consisted of several sharp and well spaced out peaks. The
intensity is not largest in the centre, but the maxima occur at around y=(+/-)150mm instead, unlike for the slit.

Supplementary Task ("CoreTasks.py"):
The aperture function was modified so that the calculation is accurate in the near-field by multiplying it by a phase
correction term. The previous calculations were repeated at large screen distances D, and found to agree with previous
Fraunhofer calculations, as expected. Then, the diffraction patterns were calculated in the near field for both the slit
and the phase grating. The intensity pattern for the slit aperture is similar to before, but not a pure sinc^2 function.
This is to be expected due to near-field phase variations across the aperture. The intensity pattern for the phase
grating is also broadly similar to the Fraunhofer pattern, with the same number of peaks and relative intensity, however
the signal is now asymmetric, which is non-physical. This is likely a result of the discrete sampling of the aperture
function. Increasing the number of samples did not improve the situation.
