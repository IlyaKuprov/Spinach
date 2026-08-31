# examples/quantum_tech/diamond_defects/diamond_gev0_epr_xw.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/diamond_defects/diamond_gev0_epr_xw.m`
- Signature: `diamond_gev0_epr_xw()`
- Total lines: 67

## Purpose

Field-swept powder EPR spectra of GeV0 centre in diamond at X and W bands. Calculation time: seconds.

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Set GeV0 centre model parameters.; implemented by `gev0_params.orientation='111'`.
- Lines 14-15: Build the spin system.; implemented by `[sys,inter]=diamond_gev0(gev0_params)`.
- Lines 17-18: Field sweep; implemented by `sys.magnet=1`.
- Lines 20-21: Define the basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 24-25: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: Set common EPR parameters; implemented by `parameters.spins={'E3'}`.
- Lines 37-38: Set X-band parameters; implemented by `parameters.mw_freq=9.5e9`.
- Lines 41-42: Run the X-band simulation; implemented by `[spec_x,par_x]=fieldsweep(spin_system,parameters)`.
- Lines 44-45: Plot the X-band spectrum; implemented by `kfigure(); scale_figure([1.50 0.75])`.
- Lines 52-53: Set W-band parameters; implemented by `parameters.mw_freq=94e9`.
- Lines 56-57: Run the W-band simulation; implemented by `[spec_w,par_w]=fieldsweep(spin_system,parameters)`.
- Lines 59-60: Plot the W-band spectrum; implemented by `subplot(1,2,2); plot(par_w.b_axis,spec_w)`.

### Key state/data transformations

- Lines 11: computes `gev0_params.orientation` using `gev0_params.orientation='111'`.
- Lines 12: computes `gev0_params.germanium` using `gev0_params.germanium='none'`.
- Lines 15: computes `[sys,inter]` using `[sys,inter]=diamond_gev0(gev0_params)`.
- Lines 18: computes `sys.magnet` using `sys.magnet=1`.
- Lines 21: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 29: computes `parameters.spins` using `parameters.spins={'E3'}`.
- Lines 30: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 31: computes `parameters.fwhm` using `parameters.fwhm=1e-3`.
- Lines 32: computes `parameters.int_tol` using `parameters.int_tol=1e-3`.
- Lines 33: computes `parameters.tm_tol` using `parameters.tm_tol=0.01`.
- Lines 34: computes `parameters.npoints` using `parameters.npoints=2048`.
- Lines 35: computes `parameters.rspt_order` using `parameters.rspt_order=Inf`.
- Lines 38: computes `parameters.mw_freq` using `parameters.mw_freq=9.5e9`.
- Lines 39: computes `parameters.window` using `parameters.window=[0.05 0.45]`.
- Lines 42: computes `[spec_x,par_x]` using `[spec_x,par_x]=fieldsweep(spin_system,parameters)`.
- Lines 57: computes `[spec_w,par_w]` using `[spec_w,par_w]=fieldsweep(spin_system,parameters)`.

## Implementation structure

- Field-swept powder EPR spectra of GeV0 centre
- in diamond at X and W bands.
- Calculation time: seconds.
- Set GeV0 centre model parameters.
- Build the spin system.
- Field sweep
- Define the basis set
- Run Spinach housekeeping
- Set common EPR parameters
- Set X-band parameters
- Run the X-band simulation
- Plot the X-band spectrum

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `diamond_gev0()`, `create()`, `basis()`, `fieldsweep()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`.
