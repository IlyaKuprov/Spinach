# examples/parahydrogen/altadena_propanal.m

- Signature: `altadena_propanal()`

## Purpose

Simulates the ALTADENA experiment for the parahydrogenation of acrolein into propanal. The simple ALTADENA model assumes perfectly adiabatic transfer and ignores isotropic mixing at low field. Note the small flip angle. Calculation time: seconds.

## Physical / mathematical content

- The six-spin `1H` system uses a magnetic field of 7.05, chemical shifts of 1.11, 1.11, 1.11, 2.46, 2.46, and 9.79, and scalar couplings of 7.3 between each of spins 1–3 and spins 4–5, and 1.4 between spin 6 and each of spins 4–5.
- The initial state is `1.0*state(spin_system,{'Lz','Lz'},{1,4}) - 0.5*state(spin_system,{'Lz'},{1}) + 0.5*state(spin_system,{'Lz'},{4})`. Detection uses `state(spin_system,'L+','1H')`; the pulse operator is `operator(spin_system,'Ly','1H')` with a small flip angle of `pi/100`.

## Numerical / algorithmic content

- Uses the `sphten-liouv` formalism with no basis approximation and symmetry groups `S3` on spins `[1 2 3]` and `S2` on spins `[4 5]`.
- Acquires the signal with `liquid(spin_system,@hp_acquire,parameters,'nmr')`, applies exponential apodisation with parameter 6, and computes `fftshift(fft(fid,parameters.zerofill))`.
- Acquisition and plotting parameters are `decouple={}`, `offset=500`, `sweep=1000`, `npoints=1024`, `zerofill=8192`, `axis_units='ppm'`, and `invert_axis=1`. The real spectrum is plotted with `plot_1d`.

## Implementation structure

- Defines the spin system, magnetic field, chemical shifts, scalar couplings, and basis; creates and configures the Spinach spin system; sets the sequence parameters; then acquires, apodises, Fourier-transforms, and plots the signal.

Source authors: Ronghui Zhou (hui@ufl.edu) and Ilya Kuprov (ilya.kuprov@weizmann.ac.il).