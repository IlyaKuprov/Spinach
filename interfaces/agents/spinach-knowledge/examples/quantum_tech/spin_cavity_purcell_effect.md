# examples/quantum_tech/spin_cavity_purcell_effect.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/spin_cavity_purcell_effect.m`
- Signature: `spin_cavity_purcell_effect()`
- Total lines: 124

## Purpose

Cavity-induced spin relaxation in the EPR Purcell regime. Coherent Jaynes-Cummings exchange is combined with rapid cavity damping in Liouville space, producing relaxation of the spin excitation by the NMR mechanism known as relaxation of the second kind. Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The physics is Jaynes-Cummings-like cavity QED: a two-level or few-level matter degree of freedom exchanges excitations with a quantised harmonic mode through rotating terms such as a†σ_- + aσ_+.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file also defines local helper function(s): `purcell_device()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Magnet field; implemented by `sys.magnet=0`.
- Lines 16-17: Particle specification; implemented by `sys.isotopes={'E','C3'}`.
- Lines 19-20: Formalism and basis; implemented by `bas.formalism='zeeman-liouv'`.
- Lines 23-24: Purcell parameters; implemented by `coupling=0.35e6`.
- Lines 28-29: Preallocate rate array; implemented by `rates=zeros(numel(detuning),numel(loss_rates))`.
- Lines 31-32: Loop over cavity loss rates; implemented by `for n=1:numel(loss_rates)`.
- Lines 34-35: Generators of the resonant damped spin-cavity device; implemented by `[spin_system,Hjc,R]=purcell_device(sys,bas,coupling,loss_rates(n))`.
- Lines 37-38: Spin detuning superoperator; implemented by `spin_ham=operator(spin_system,'Lz',1)`.
- Lines 40-41: Extract Purcell rates from Liouvillian slow modes; implemented by `for k=1:numel(detuning)`.
- Lines 43-44: Assemble the dissipative Liouvillian at this detuning; implemented by `L=-1i*(Hjc+detuning(k)*spin_ham)+R`.
- Lines 46-47: Convert the slow spin-amplitude mode into a population rate; implemented by `decay_modes=eig(full(L))`.
- Lines 55-56: Validate resonant second-kind relaxation; implemented by `if rates(ceil(end/2),2)<=rates(1,2)`.
- Lines 60-61: Validate the resonant rate against the exact two-level Purcell rate; implemented by `g_rad=2*pi*coupling; kappa_ref=loss_rates(2)`.
- Lines 67-68: Rebuild the generators at the reference loss rate; implemented by `[spin_system,Hjc,R]=purcell_device(sys,bas,coupling,loss_rates(2))`.
- Lines 71-72: Pick representative detunings for survival curves; implemented by `time_axis=linspace(0,40e-6,250)`.
- Lines 78-79: Simulate spin excitation survival under cavity damping; implemented by `for n=1:numel(det_pick)`.
- Lines 81-82: Assemble the dissipative Liouvillian at this detuning; implemented by `L=-1i*(Hjc+det_pick(n)*spin_ham)+R`.
- Lines 84-85: Propagate the density matrix in Liouville space; implemented by `for k=1:numel(time_axis)`.

### Control flow inferred from the code

- Line 32: `for` loop over `n=1:numel(loss_rates)`.
- Line 41: `for` loop over `k=1:numel(detuning)`.
- Line 56: conditional branch on `rates(ceil(end/2),2)<=rates(1,2)`.
- Line 63: conditional branch on `abs(rates(ceil(end/2),2)-purcell)>1e-4*purcell`.
- Line 79: `for` loop over `n=1:numel(det_pick)`.
- Line 85: `for` loop over `k=1:numel(time_axis)`.
- Line 92: conditional branch on `survival(end,1)>=survival(end,end)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=0`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'E','C3'}`.
- Lines 20: computes `bas.formalism` using `bas.formalism='zeeman-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `coupling` using `coupling=0.35e6`.
- Lines 25: computes `detuning` using `detuning=2*pi*linspace(-8e6,8e6,301)`.
- Lines 26: computes `loss_rates` using `loss_rates=2*pi*[2e6 4e6 8e6 16e6]`.
- Lines 29: computes `rates` using `rates=zeros(numel(detuning),numel(loss_rates))`.
- Lines 35: computes `[spin_system,Hjc,R]` using `[spin_system,Hjc,R]=purcell_device(sys,bas,coupling,loss_rates(n))`.
- Lines 38: computes `spin_ham` using `spin_ham=operator(spin_system,'Lz',1)`.
- Lines 44: computes `L` using `L=-1i*(Hjc+detuning(k)*spin_ham)+R`.
- Lines 47: computes `decay_modes` using `decay_modes=eig(full(L))`.
- Lines 49: computes `rates(k,n)` using `rates(k,n)=-2*max(real(decay_modes))`.
- Lines 61: computes `g_rad` using `g_rad=2*pi*coupling; kappa_ref=loss_rates(2)`.
- Lines 62: computes `purcell` using `purcell=kappa_ref/2-sqrt(kappa_ref^2/4-4*g_rad^2)`.
- Lines 72: computes `time_axis` using `time_axis=linspace(0,40e-6,250)`.
- Lines 73: computes `det_pick` using `det_pick=2*pi*[0 2e6 6e6]`.
- Lines 74: computes `survival` using `survival=zeros(numel(time_axis),numel(det_pick))`.

### Local helper functions

- Line 113: `purcell_device()` — `function [spin_system,Hjc,R]=purcell_device(sys,bas,coupling,loss_rate)`.
  - Representative operation: `inter.modes.frqs={[] 0}`.
  - Representative operation: `inter.modes.linewidths={[] loss_rate/(2*pi)}`.

## Implementation structure

- Cavity-induced spin relaxation in the EPR Purcell regime.
- Coherent Jaynes-Cummings exchange is combined with rapid
- cavity damping in Liouville space, producing relaxation of
- the spin excitation by the NMR mechanism known as relaxation
- of the second kind.
- Calculation time: seconds
- Magnet field
- Particle specification
- Formalism and basis
- Purcell parameters
- Preallocate rate array
- Loop over cavity loss rates

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `purcell_device()`, `loss_rates()`, `operator()`, `detuning()`, `decay_modes()`, `rates()`, `state()`, `det_pick()`, `survival()`, `time_axis()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`.
