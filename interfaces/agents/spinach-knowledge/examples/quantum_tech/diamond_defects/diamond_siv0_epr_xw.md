# examples/quantum_tech/diamond_defects/diamond_siv0_epr_xw.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/diamond_defects/diamond_siv0_epr_xw.m`
- Signature: `diamond_siv0_epr_xw()`
- Total lines: 68

## Purpose

Field-swept powder EPR spectra of SiV0 centre in diamond at X and W bands. Calculation time: seconds.

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Set SiV0 centre model parameters.; implemented by `siv0_params.orientation='111'`.
- Lines 15-16: Build the spin system; implemented by `[sys,inter]=diamond_siv0(siv0_params)`.
- Lines 18-19: Field sweep; implemented by `sys.magnet=1`.
- Lines 21-22: Define the basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 25-26: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 29-30: EPR sim parameters; implemented by `parameters.spins={'E3'}`.
- Lines 38-39: Set X-band parameters; implemented by `parameters.mw_freq=9.5e9`.
- Lines 42-43: Run the X-band simulation; implemented by `[spec_x,par_x]=fieldsweep(spin_system,parameters)`.
- Lines 45-46: Plot the X-band spectrum; implemented by `kfigure(); scale_figure([1.50 0.75])`.
- Lines 53-54: Set W-band parameters; implemented by `parameters.mw_freq=94e9`.
- Lines 57-58: Run the W-band simulation; implemented by `[spec_w,par_w]=fieldsweep(spin_system,parameters)`.
- Lines 60-61: Plot the W-band spectrum; implemented by `subplot(1,2,2); plot(par_w.b_axis,spec_w)`.

### Key state/data transformations

- Lines 11: computes `siv0_params.orientation` using `siv0_params.orientation='111'`.
- Lines 12: computes `siv0_params.silicon` using `siv0_params.silicon='29Si'`.
- Lines 13: computes `siv0_params.n_13c` using `siv0_params.n_13c=0`.
- Lines 16: computes `[sys,inter]` using `[sys,inter]=diamond_siv0(siv0_params)`.
- Lines 19: computes `sys.magnet` using `sys.magnet=1`.
- Lines 22: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 30: computes `parameters.spins` using `parameters.spins={'E3'}`.
- Lines 31: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 32: computes `parameters.fwhm` using `parameters.fwhm=0.001`.
- Lines 33: computes `parameters.int_tol` using `parameters.int_tol=1e-4`.
- Lines 34: computes `parameters.tm_tol` using `parameters.tm_tol=0.001`.
- Lines 35: computes `parameters.npoints` using `parameters.npoints=2048`.
- Lines 36: computes `parameters.rspt_order` using `parameters.rspt_order=Inf`.
- Lines 39: computes `parameters.mw_freq` using `parameters.mw_freq=9.5e9`.
- Lines 40: computes `parameters.window` using `parameters.window=[0.1 0.5]`.
- Lines 43: computes `[spec_x,par_x]` using `[spec_x,par_x]=fieldsweep(spin_system,parameters)`.

## Implementation structure

- Field-swept powder EPR spectra of SiV0 centre
- in diamond at X and W bands.
- Calculation time: seconds.
- Set SiV0 centre model parameters.
- Build the spin system
- Field sweep
- Define the basis set
- Run Spinach housekeeping
- EPR sim parameters
- Set X-band parameters
- Run the X-band simulation
- Plot the X-band spectrum

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `diamond_siv0()`, `create()`, `basis()`, `fieldsweep()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`.
