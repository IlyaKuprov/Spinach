# examples/relaxation_theory/csa_csa_xcorr_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/csa_csa_xcorr_2.m`
- Signature: `csa_csa_xcorr_2()`
- Total lines: 63

## Purpose

CSA-CSA cross-correlation in the 103Rh subsystem and its effect on the widths of the three lines of the proton triplet. Calculation time: seconds.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnet field; implemented by `sys.magnet=11.75`.
- Lines 13-14: Set the spin system; implemented by `sys.isotopes={'1H','103Rh','103Rh'}`.
- Lines 19-20: J-couplings; implemented by `inter.coupling.scalar={0 4 4`.
- Lines 24-25: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Sequence parameters -1H; implemented by `parameters.spins={'1H'}`.
- Lines 50-51: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 53-54: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',20}})`.
- Lines 56-57: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 59-60: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=11.75`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','103Rh','103Rh'}`.
- Lines 15-17: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={[6.9 0 0; 0 6.9 0; 0 0 6.9] [7250 0 0; 0 8000 0; 0 0 7250] [7250 0 0; 0 8000 0; 0 0 7250]}`.
- Lines 20: computes `inter.coupling.scalar` using `inter.coupling.scalar={0 4 4`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 27: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 28: computes `inter.tau_c` using `inter.tau_c={10e-9}`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 40: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 41: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 42: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 43: computes `parameters.offset` using `parameters.offset=6.9*500`.
- Lines 44: computes `parameters.sweep` using `parameters.sweep=50`.
- Lines 45: computes `parameters.npoints` using `parameters.npoints=2048`.

## Implementation structure

- CSA-CSA cross-correlation in the 103Rh subsystem and its effect
- on the widths of the three lines of the proton triplet.
- Calculation time: seconds.
- Magnet field
- Set the spin system
- J-couplings
- Relaxation theory
- Basis set
- Spinach housekeeping
- Sequence parameters -1H
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
