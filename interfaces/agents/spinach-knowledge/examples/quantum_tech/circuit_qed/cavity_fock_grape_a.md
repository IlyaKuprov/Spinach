# examples/quantum_tech/circuit_qed/cavity_fock_grape_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/circuit_qed/cavity_fock_grape_a.m`
- Signature: `cavity_fock_grape_a()`
- Total lines: 127

## Purpose

GRAPE preparation of a cavity Fock state through a dispersively coupled qubit, using piecewise-constant drives on both the cavity and the qubit. A linear drive alone cannot make a Fock state out of the vacuum of a harmonic mode; the qubit conditions the cavity phase through the dispersive shift and thereby provides the requi- red nonlinearity. The optimisation is run at two Fock space trun- cations; the pulse optimis

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file also defines local helper function(s): `pulse_fidelity()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Cavity truncation levels to compare; implemented by `trunc_labels={'C3','C4'}`.
- Lines 22-23: Preallocate pulse and fidelity storage; implemented by `pulses=cell(1,2); fids=zeros(1,3)`.
- Lines 25-26: Loop over the cavity truncations; implemented by `for m=1:2`.
- Lines 28-29: Magnet field; implemented by `sys.magnet=0`.
- Lines 31-32: Truncated cavity mode and a qubit; implemented by `sys.isotopes={trunc_labels{m},'E'}`.
- Lines 34-35: Cavity on resonance with its drive; implemented by `inter.modes.frqs={0 []}`.
- Lines 37-38: Dispersive coupling to the qubit; implemented by `inter.modes.dispersive=cell(2,2)`.
- Lines 41-42: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 45-46: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 49-50: Drift Hamiltonian from the declared interactions; implemented by `H=hamiltonian(assume(spin_system,'cavity'))`.
- Lines 52-53: Quadrature control operators of the cavity and the qubit; implemented by `Cr=operator(spin_system,'C',1); An=operator(spin_system,'A',1)`.
- Lines 57-58: Vacuum cavity and cavity Fock state two, qubit in the upper level; implemented by `rho_init=state(spin_system,{'BL1','ZL2'},{1,2})`.
- Lines 63-64: Define control parameters; implemented by `control.isotopes={trunc_labels{m},'E'}`.
- Lines 77-78: Plots during optimisation; implemented by `control.plotting={'xy_controls'}`.
- Lines 80-81: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 83-84: Gaussian initial guess on the in-phase channels; implemented by `env=exp(-(33e-9*((1:40)-0.5)-0.66e-6).^2/(2*(1.32e-6/8)^2))`.
- Lines 87-88: Run the optimisation, get normalised pulse; implemented by `pulses{m}=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 90-91: Recompute the transfer fidelity by direct propagation; implemented by `fids(m)=pulse_fidelity(spin_system,H,ops,pulses{m},rho_init,rho_targ)`.

### Control flow inferred from the code

- Line 26: `for` loop over `m=1:2`.
- Line 99: conditional branch on `(fids(1)<0.95)||(fids(2)<0.95)`.

### Key state/data transformations

- Lines 20: computes `trunc_labels` using `trunc_labels={'C3','C4'}`.
- Lines 23: computes `pulses` using `pulses=cell(1,2); fids=zeros(1,3)`.
- Lines 29: computes `sys.magnet` using `sys.magnet=0`.
- Lines 32: computes `sys.isotopes` using `sys.isotopes={trunc_labels{m},'E'}`.
- Lines 35: computes `inter.modes.frqs` using `inter.modes.frqs={0 []}`.
- Lines 38: computes `inter.modes.dispersive` using `inter.modes.dispersive=cell(2,2)`.
- Lines 39: computes `inter.modes.dispersive{1,2}` using `inter.modes.dispersive{1,2}=656.2e3`.
- Lines 42: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 43: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 46: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 50: computes `H` using `H=hamiltonian(assume(spin_system,'cavity'))`.
- Lines 53: computes `Cr` using `Cr=operator(spin_system,'C',1); An=operator(spin_system,'A',1)`.
- Lines 54: computes `Lx` using `Lx=operator(spin_system,'Lx',2); Ly=operator(spin_system,'Ly',2)`.
- Lines 55: computes `ops` using `ops={(Cr+An)/2,1i*(Cr-An)/2,Lx,Ly}`.
- Lines 58: computes `rho_init` using `rho_init=state(spin_system,{'BL1','ZL2'},{1,2})`.
- Lines 59: computes `rho_targ` using `rho_targ=state(spin_system,{'BL3','ZL2'},{1,2})`.
- Lines 64: computes `control.isotopes` using `control.isotopes={trunc_labels{m},'E'}`.
- Lines 65: computes `control.channels` using `control.channels=[1;1;2;2]`.

### Local helper functions

- Line 116: `pulse_fidelity()` — `function fid=pulse_fidelity(spin_system,H,ops,pulse,rho_init,rho_targ)`.
  - Representative operation: `rho=rho_init`.
  - Representative operation: `for n=1:size(pulse,2)`.

## Implementation structure

- GRAPE preparation of a cavity Fock state through a dispersively
- coupled qubit, using piecewise-constant drives on both the cavity
- and the qubit. A linear drive alone cannot make a Fock state out
- of the vacuum of a harmonic mode; the qubit conditions the cavity
- phase through the dispersive shift and thereby provides the requi-
- red nonlinearity. The optimisation is run at two Fock space trun-
- cations; the pulse optimised in the smaller space underperforms
- when it is re-evaluated in the larger one -optimal control solu-
- tions must be converged with respect to the Fock space truncation.
- Model and parameters from the bosonic GRAPE example of the para-
- qeet package.
- Calculation time: minutes

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `state()`, `optimcon()`, `fmaxnewton()`, `fids()`, `pulse_fidelity()`, `num2str()`, `kfigure()`, `bar()`, `set()`, `kylabel()`, `ktitle()`.
