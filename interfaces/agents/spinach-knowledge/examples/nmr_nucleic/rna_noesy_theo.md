# examples/nmr_nucleic/rna_noesy_theo.m

- Signature: `rna_noesy_theo()`

## Purpose

Simulates a 1H–1H NOESY spectrum of the example RNA provided by Gerhard Wagner's group at Harvard University. The source estimates a calculation time of hours and credits Shunsuke Imai, Scott Robson, Gerhard Wagner, Zenawi Welderufael, and Ilya Kuprov.

## Physical / mathematical content

The RNA structure and assignments are imported from `example.pdb` and `example.txt`; the listed exchangeable protons are deuterated, and 13C and 15N spins are removed under the source's assumption that the RNA is unlabelled. The system is set to 17.62 T with Redfield relaxation (`rlx_keep='kite'`, zero equilibrium, and `tau_c={3e-9}`). The basis uses the sphten-liouv formalism, IK-1 approximation, scalar couplings, interaction level 5, and proximity level 3.

## Numerical / algorithmic content

The sequence uses a 0.200 mixing time, 1H spins, offsets of 3473, sweeps [7500, 7500], 512 acquired points, and zero filling [1024, 4096]. The source disables Krylov propagation and enables propagator caching and the greedy option. It applies squared-cosine apodisation, Fourier-transforms F2, combines the cosine and sine components as a States signal, Fourier-transforms F1, then plots the negative real spectrum.

## Implementation structure

After system construction, the code removes 13C and 15N spins, builds the basis, and simulates with `liquid(spin_system,@noesy,parameters,'nmr')`. The processed spectrum is displayed with the source's plotting scale and contour settings.
