# experiments/slowpass_sat.m

Computes a saturated continuous-wave spectrum by solving the stationary driven Liouvillian at each offset. The irradiation is retained in the generator rather than treated as a linear-response perturbation. Use the five-argument context callback signature shared with slowpass, but not its amplitude or preparation conventions.

parameters.irr_powers is a finite non-negative scalar nutation frequency nu_1 in Hz; the drive is 2*pi*nu_1*Lx. parameters.spins selects one isotope, parameters.sweep contains two finite offsets in Hz, parameters.npoints is an integer greater than one, and parameters.coil is one finite detection column. The returned complex column is the detector expectation, not an FFT-scaled spectrum or absorbed power. rho0 is ignored. A detector containing an identity component includes the known pinned contribution; use coil_state for concentration-independent detection.

Requires sphten-liouv, explicit matrices, thermalised relaxation, and a unique stationary solution after each substance population is fixed. Population-changing kinetics are rejected. Spin-only contexts and gridfree with its uniform orientational distribution are supported; other spatial contexts and orienting potentials are outside the contract. Keep longitudinal and both transverse coherence orders in the basis. No trajectory reduction is applied.

The CPU path uses direct backslash up to 100000 active coordinates and scaled, shifted-ILU right-preconditioned GMRES above that. Factors are recomputed at each offset; the shift affects only the preconditioner. Iterative convergence and the unpreconditioned relative residual are checked at 1e-6. There is no GPU or polyadic path.

The three sle_esr_sat examples in examples/esr_sol_swept illustrate nitroxide microwave saturation and rotational diffusion using the parameters of Gupta, Dzikovski, and Freed, DOI 10.1007/s00723-021-01390-7. Their headers distinguish the adopted relaxation, sampling, and amplitude conventions from exact published numerical reproduction.
