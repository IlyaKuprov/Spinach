# experiments/spen/spendosycosy.m

[Canonical source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/spen/spendosycosy.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=spendosycosy.m)

Ultrafast 3D DOSY-COSY combines the DOSY preparation and diffusion delay with a COSY evolution dimension. The source forms L = H + F + 1i*R + 1i*K, excites rho0 with pi/2, and steps through the coherence orders +1 → -1 → 0 during the first chirp/encoding interval and diffusion preparation. After the delay td-Tau-Te, a further pi/2 pulse takes the state through -1; a second chirp and positive Ge*G{1} interval return it to +1. The chirp parameters are pulsenpoints, Te, BW, smfactor, and chirptype ('wurst' or 'smoothed').

The COSY axis is generated as a trajectory under L with timestep 1/sweep and npoints2-1 intervals, then a pi/2 pulse and +1 selection are applied to that stack. The acquisition trace uses negative-gradient prephasing for half of npoints1*deltat, followed by npoints1 detections of coil'*rho under the positive acquisition gradient. The returned fid is [npoints1, npoints2, nloops]: acquisition points, COSY evolution points, and loop index. npoints2 is the trajectory length, including its initial state. The header also mentions singular `npoints`, but the body reads and validates `npoints1` and `npoints2` instead. Loop bodies run with parfor; the source transfers the COSY stack, propagators, and coil to the GPU when GPU execution is enabled.

The function checks rho0, coil, scalar dims (m), scalar npts, spins, npoints1, npoints2, sweep, deltat, nloops, Ga, pulsenpoints, smfactor, Te, Tau, td, BW, Ge, and chirptype in parameters. Gradient amplitudes are in T/m, and td must be at least Tau+Te. Although the source header lists Gp and Tp as coherence-selection gradient settings, the function body does not access those fields. The required formalism is sphten-liouv; H, R, K, and F must be equal-sized matrices, G must be a cell array, and these context inputs are supplied by imaging.

Authors: jeannicolas.dumez@cnrs.fr, ilya.kuprov@weizmann.ac.il, ludmilla.guduff@cnrs.fr.
