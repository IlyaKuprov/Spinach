# experiments/spen/spencosy.m

[Canonical source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/spen/spencosy.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=spencosy.m)

Ultrafast COSY with spatial encoding and gradient acquisition. The source forms the Fokker–Planck generator L = H + F + 1i*R + 1i*K and applies an initial pi/2 pulse to parameters.rho0. It then applies two shaped chirp pulses from chirp_pulse with opposite signs of the encoding gradient Ge*G{1}. A gradient interval Gp*G{1} for Tp brackets the next pi/2 pulse. The source does not insert explicit coherence-order projection calls in this COSY sequence.

The trace starts after negative acquisition-gradient prephasing for `Taq/2`, where `Taq=parameters.npoints*parameters.deltat`. Each loop-start state is then advanced by one combined positive/negative Ga*G{1} loop propagator; the acquired trace samples coil'*rho while stepping under the positive acquisition gradient. fid has shape [parameters.npoints, parameters.nloops]: points along each readout, then loop index. The loop bodies run with parfor; when GPU execution is enabled, the source moves propagators, state, and coil to the GPU.

Required sequence settings checked by the source are parameters.rho0, parameters.coil, scalar parameters.dims (sample size in m), parameters.npts (spin-packet count), parameters.spins, parameters.deltat, parameters.npoints, parameters.nloops, parameters.Ga, parameters.pulsenpoints, parameters.nWURST, parameters.Te, parameters.BW, parameters.Ge, parameters.Gp, and parameters.Tp. The source header gives gradient amplitudes in T/m and identifies nWURST as a pulse-smoothing factor. The header also documents D in m^2/s, but the function body does not read parameters.D; diffusion/flow enters through F. H, R, K, and F must be equal-sized matrices, G must be a cell array, and the required formalism is sphten-liouv. H, R, K, G, and F are supplied by the imaging context.

Authors: jeannicolas.dumez@cnrs.fr, ilya.kuprov@weizmann.ac.il, ludmilla.guduff@cnrs.fr.
