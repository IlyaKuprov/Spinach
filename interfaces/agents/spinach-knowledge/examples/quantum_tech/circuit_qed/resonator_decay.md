# examples/quantum_tech/circuit_qed/resonator_decay.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/circuit_qed/resonator_decay.m`
- Signature: `resonator_decay()`
- Total lines: 117

## Purpose

Open-system dynamics of a leaky microwave resonator at finite temperature. A Fock state decays as a downward cascade through the level ladder, and a coherent state decays with its Poisson population structure largely preserved; both settle into the thermal state of the mode. The mean photon number follows the analytical amplitude damping solution in both cases, up to the distortion of the weak thermal channel by the 

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Magnet field; implemented by `sys.magnet=0`.
- Lines 20-21: Microwave resonator with five Fock levels; implemented by `sys.isotopes={'C5'}`.
- Lines 23-24: Resonator frequency; implemented by `inter.modes.frqs={6.02e9}`.
- Lines 26-27: Energy relaxation lifetime; implemented by `inter.modes.lifetimes={10e-9}`.
- Lines 29-30: Temperature of the environment; implemented by `inter.temperature=0.050`.
- Lines 32-34: Bose-Einstein occupation at the environment temperature; implemented by `n_eq=1/(exp(6.62607015e-34*inter.modes.frqs{1}/ (1.380649e-23*inter.temperature))-1)`.
- Lines 36-37: Total coherence time from a five microsecond pure dephasing time; implemented by `inter.modes.t2_times={1/(1/5e-6+(1+2*n_eq)/(2*inter.modes.lifetimes{1}))}`.
- Lines 39-40: Formalism and basis; implemented by `bas.formalism='zeeman-liouv'`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Dissipative evolution generator; implemented by `H=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 52-53: Fock state and coherent state as initial conditions; implemented by `rho_fock=state(spin_system,'BL5',1)`.
- Lines 56-59: Fock level population detection states; implemented by `coils=[state(spin_system,'BL1',1) state(spin_system,'BL2',1) state(spin_system,'BL3',1) state(spin_system,'BL4',1) state(spin_system,'BL5',1)]`.
- Lines 61-62: Time grid of the source calculation; implemented by `nsteps=100; dt=1e-9; time_axis=dt*(0:nsteps)`.
- Lines 64-65: Propagator over one time step; implemented by `P=expm(full(G)*dt)`.
- Lines 67-68: Preallocate population trajectories; implemented by `pops_fock=zeros(5,nsteps+1); pops_coh=zeros(5,nsteps+1)`.
- Lines 70-71: Propagate and record Fock level populations; implemented by `for n=1:(nsteps+1)`.
- Lines 77-78: Mean photon numbers from the populations; implemented by `n_fock=(0:4)*pops_fock; n_coh=(0:4)*pops_coh`.
- Lines 80-81: Validate against the analytical solution, Fock truncation limits accuracy; implemented by `kappa=1/inter.modes.lifetimes{1}`.

### Control flow inferred from the code

- Line 71: `for` loop over `n=1:(nsteps+1)`.
- Line 84: conditional branch on `max(abs(n_fock-n_fock_ref))>5e-3`.
- Line 87: conditional branch on `max(abs(n_coh-n_coh_ref))>5e-3`.
- Line 94: conditional branch on `~((idx_three<idx_two)&&(idx_two<idx_one))`.
- Line 99: conditional branch on `abs(pops_fock(1,end)-1/(1+n_eq))>2e-3`.

### Key state/data transformations

- Lines 18: computes `sys.magnet` using `sys.magnet=0`.
- Lines 21: computes `sys.isotopes` using `sys.isotopes={'C5'}`.
- Lines 24: computes `inter.modes.frqs` using `inter.modes.frqs={6.02e9}`.
- Lines 27: computes `inter.modes.lifetimes` using `inter.modes.lifetimes={10e-9}`.
- Lines 30: computes `inter.temperature` using `inter.temperature=0.050`.
- Lines 33-34: computes `n_eq` using `n_eq=1/(exp(6.62607015e-34*inter.modes.frqs{1}/ (1.380649e-23*inter.temperature))-1)`.
- Lines 37: computes `inter.modes.t2_times` using `inter.modes.t2_times={1/(1/5e-6+(1+2*n_eq)/(2*inter.modes.lifetimes{1}))}`.
- Lines 40: computes `bas.formalism` using `bas.formalism='zeeman-liouv'`.
- Lines 41: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 44: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 48: computes `H` using `H=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 49: computes `R` using `R=relaxation(spin_system)`.
- Lines 50: computes `G` using `G=-1i*H+R`.
- Lines 53: computes `rho_fock` using `rho_fock=state(spin_system,'BL5',1)`.
- Lines 54: computes `rho_coh` using `rho_coh=coherent(spin_system,1,1.5)`.
- Lines 57-59: computes `coils` using `coils=[state(spin_system,'BL1',1) state(spin_system,'BL2',1) state(spin_system,'BL3',1) state(spin_system,'BL4',1) state(spin_system,'BL5',1)]`.
- Lines 62: computes `nsteps` using `nsteps=100; dt=1e-9; time_axis=dt*(0:nsteps)`.
- Lines 65: computes `P` using `P=expm(full(G)*dt)`.

## Implementation structure

- Open-system dynamics of a leaky microwave resonator at finite
- temperature. A Fock state decays as a downward cascade through
- the level ladder, and a coherent state decays with its Poisson
- population structure largely preserved; both settle into the
- thermal state of the mode. The mean photon number follows the
- analytical amplitude damping solution in both cases, up to the
- distortion of the weak thermal channel by the Fock space trun-
- cation. Model and parameters from the resonator decay example
- of the paraqeet package.
- Calculation time: seconds
- Magnet field
- Microwave resonator with five Fock levels

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `relaxation()`, `state()`, `coherent()`, `pops_fock()`, `pops_coh()`, `n_fock()`, `n_coh()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`.
