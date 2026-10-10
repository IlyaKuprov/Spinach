# examples/relaxation_theory/cpmg_echo_train.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/cpmg_echo_train.m)

This script computes a powder-averaged CPMG echo train for two protons at 14.1 T, then plots its real time-domain signal. It is a model calculation, not a fit to or reproduction of a named experiment. The two shielding principal-value expressions are [-2 -2 4]-5 and [-1 -3 4]+5 (ppm), and both Euler-angle rows are zero. The example uses the phenomenological t1_t2 relaxation model with per-spin R1 rates {50.0 50.0} and R2 rates {150.0 150.0}; the source does not annotate units for these rate inputs. Equilibrium is set to zero and relaxation retention is secular. The basis is sphten-liouv without approximation, and trajlevel is disabled.

The powder grid is rep_2ang_200pts_sph. The selected channel is 1H, with L+ for both the initial state and receiver coil and Lx as the pulse operator. The parameters set 10 loops, a 1e-5 s timestep and 100 points. The powder NMR simulation returns fid; the plotted observable is real(fid) against the constructed time axis, labelled as the S_X expectation value. No measured echo train is supplied for comparison.
