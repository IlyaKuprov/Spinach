# examples/esr_sol_swept/sle_esr_sat_rate.m

Demonstrates the rotational-rate dependence of saturated X-band nitroxide spectra at B1=0.5 G. The rates 10^[7,7.5,8,9] per second and magnetic tensors follow Figure 3 of Gupta, Dzikovski, and Freed, DOI 10.1007/s00723-021-01390-7. The isotropic second-rank correlation time is 1/(6R).

Uses electron R1=2e5 and R2=1e5 per second and IME thermalisation at 298 K. Plots the field derivative of minus the imaginary transverse expectation over 3300–3420 G, centred at 3360 G through the isotropic-g frequency reference. The panels retain their relative thermal-magnetisation amplitude scale, not independent peak normalisation. Basis truncation and field sampling must be converged for quantitative use; the example does not include an orienting potential.
