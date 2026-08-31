# examples/quantum_tech/diamond_defects/diamond_p1_13c_epr_xw.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/diamond_defects/diamond_p1_13c_epr_xw.m`
- Signature: `diamond_p1_13c_epr_xw()`
- Total lines: 71

## Purpose

Field-swept powder EPR spectra of a P1 centre in 13C-enriched diamond at X and W bands. Calculation time: minutes.

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Set P1 model parameters.; implemented by `p1_params.orientation='111'`.
- Lines 14-15: Build the spin system.; implemented by `[sys,inter]=diamond_p1_13c(p1_params)`.
- Lines 17-18: Field sweep; implemented by `sys.magnet=1`.
- Lines 20-21: Define the basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 24-25: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 27-29: Leave only 14N nucleus and 13C nuclei with hyperfines larger than a_iso > 8 MHz; implemented by `spin_system=kill_spin(spin_system,[5 7:9 11:15 17:20])`.
- Lines 33-34: Set common EPR parameters; implemented by `parameters.spins={'E'}`.
- Lines 42-43: Set X-band parameters; implemented by `parameters.mw_freq=9.5e9`.
- Lines 46-47: Run the X-band simulation; implemented by `[spec_x,par_x]=fieldsweep(spin_system,parameters)`.
- Lines 49-50: Plot the X-band spectrum; implemented by `kfigure(); scale_figure([1.50 0.75])`.
- Lines 57-58: Set W-band parameters; implemented by `parameters.mw_freq=94e9`.
- Lines 61-62: Run the W-band simulation; implemented by `[spec_w,par_w]=fieldsweep(spin_system,parameters)`.
- Lines 64-65: Plot the W-band spectrum; implemented by `subplot(1,2,2); plot(par_w.b_axis,spec_w)`.

### Key state/data transformations

- Lines 11: computes `p1_params.orientation` using `p1_params.orientation='111'`.
- Lines 12: computes `p1_params.nitrogen` using `p1_params.nitrogen='14N'`.
- Lines 15: computes `[sys,inter]` using `[sys,inter]=diamond_p1_13c(p1_params)`.
- Lines 18: computes `sys.magnet` using `sys.magnet=1`.
- Lines 21: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 35: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 36: computes `parameters.fwhm` using `parameters.fwhm=5e-4`.
- Lines 37: computes `parameters.int_tol` using `parameters.int_tol=0.1`.
- Lines 38: computes `parameters.tm_tol` using `parameters.tm_tol=0.1`.
- Lines 39: computes `parameters.npoints` using `parameters.npoints=1024`.
- Lines 40: computes `parameters.rspt_order` using `parameters.rspt_order=Inf`.
- Lines 43: computes `parameters.mw_freq` using `parameters.mw_freq=9.5e9`.
- Lines 44: computes `parameters.window` using `parameters.window=[0.31 0.36]`.
- Lines 47: computes `[spec_x,par_x]` using `[spec_x,par_x]=fieldsweep(spin_system,parameters)`.
- Lines 62: computes `[spec_w,par_w]` using `[spec_w,par_w]=fieldsweep(spin_system,parameters)`.

## Implementation structure

- Field-swept powder EPR spectra of a P1 centre
- in 13C-enriched diamond at X and W bands.
- Calculation time: minutes.
- Set P1 model parameters.
- Build the spin system.
- Field sweep
- Define the basis set
- Run Spinach housekeeping
- Leave only 14N nucleus and 13C nuclei with hyperfines larger than
- a_iso > 8 MHz
- Set common EPR parameters
- Set X-band parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `diamond_p1_13c()`, `create()`, `kill_spin()`, `basis()`, `fieldsweep()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`.
