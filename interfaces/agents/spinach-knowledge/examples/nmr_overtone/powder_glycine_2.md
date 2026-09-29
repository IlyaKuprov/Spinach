# examples/nmr_overtone/powder_glycine_2.m

- MATLAB implementation: [examples/nmr_overtone/powder_glycine_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/powder_glycine_2.m)

- Signature: `powder_glycine_2()`

## Purpose and provenance

This example calculates a powder 14N overtone NMR spectrum of glycine with Spinach's Fokker–Planck formalism. The source attributes the glycine quadrupolar tensor data to O'Dell and Ratcliffe: https://doi.org/10.1016/j.cplett.2011.08.030. That is the stated provenance of the input data; the code is a simulation, not a claim to reproduce an experiment. The source identifies T2,-2 as the coherence responsible for the overtone signal in a static sample and estimates a calculation time of seconds.

## Spin system and model

The system contains 14N at `sys.magnet=14.1`. Its quadrupolar interaction is set by `eeqq2nqi(1.18e6,0.53,1,[0 0 0])`, and the scalar Zeeman entry is `32.4`; these are the literal values and call arguments in the example. The source does not label units for these numeric inputs or further describe the conversion arguments. The basis is `sphten-liouv` with `approximation='none'`. Relaxation is `damp`, diagonal terms are retained, equilibrium is zero, and `damp_rate=500`. Krylov and trajectory-level methods are explicitly disabled in this example; a general description of Krylov-accelerated overtone propagation does not describe this setup.

## Powder acquisition and output

The source sets the magic angle to `atan(sqrt(2))`, initialises `rho0` as the 14N `T2,-2` state, and uses a coil operator formed from `cos(theta)*Lz + sin(theta)*Lx`. The powder grid is `rep_2ang_6400pts_sph`; the sweep is `[0e3 15e3]`, with 256 points and 256-point zero filling. The selected spin is 14N and the plotted axis units are kHz. The spectrum call is `powder(spin_system,@overtone_a,parameters,'qnmr')`; the plot displays its real part through `plot_1d`.

The header describes a static-sample spectrum. This script specifies no MAS or DOR rotor rate, RF pulse sequence, or contact-time parameter. No experimental comparison or fitted spectrum is produced by this file.
