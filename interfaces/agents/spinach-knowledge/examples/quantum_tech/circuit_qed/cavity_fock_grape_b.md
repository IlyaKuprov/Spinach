# examples/quantum_tech/circuit_qed/cavity_fock_grape_b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/circuit_qed/cavity_fock_grape_b.m`
- Signature: `cavity_fock_grape_b()`
- Total lines: 113

## Purpose

GRAPE preparation of a cavity Fock state through a dispersively coupled qubit using smooth band-limited drives. The controls are expanded in an orthonormal basis of slow sine and cosine waves, and the optimisation runs over the expansion coefficients, so that the resulting pulses are hardware-friendly smooth envelopes rather than free piecewise-constant switches. Compare with the piecewise-constant treatment in cavit

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Magnet field; implemented by `sys.magnet=0`.
- Lines 20-21: Truncated cavity mode and a qubit; implemented by `sys.isotopes={'C4','E'}`.
- Lines 23-24: Cavity on resonance with its drive; implemented by `inter.modes.frqs={0 []}`.
- Lines 26-27: Dispersive coupling to the qubit; implemented by `inter.modes.dispersive=cell(2,2)`.
- Lines 30-31: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Drift Hamiltonian from the declared interactions; implemented by `H=hamiltonian(assume(spin_system,'cavity'))`.
- Lines 41-42: Quadrature control operators of the cavity and the qubit; implemented by `Cr=operator(spin_system,'C',1); An=operator(spin_system,'A',1)`.
- Lines 46-47: Vacuum cavity and cavity Fock state two, qubit in the upper level; implemented by `rho_init=state(spin_system,{'BL1','ZL2'},{1,2})`.
- Lines 52-53: Band-limited orthonormal waveform basis; implemented by `wave_set=[wave_basis('sine_waves',2,40) wave_basis('cosine_waves',3,40)]'`.
- Lines 55-56: Define control parameters; implemented by `control.isotopes={'C4','E'}`.
- Lines 70-71: Plots during optimisation; implemented by `control.plotting={'xy_controls'}`.
- Lines 73-74: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 76-77: Deterministic asymmetric coefficient guess; implemented by `guess=[0.3 0.1 0.5 0.1 0.0; 0.0 0.1 0.0 0.1 0.0`.
- Lines 80-81: Run the optimisation, get basis coefficients; implemented by `coeffs=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 83-84: Reassemble the smooth time-domain pulse; implemented by `pulse=coeffs*wave_set`.
- Lines 86-87: Recompute the transfer fidelity by direct propagation; implemented by `rho=rho_init`.
- Lines 97-98: Validate the optimisation; implemented by `if fid<0.80`.

### Control flow inferred from the code

- Line 88: `for` loop over `n=1:numel(spin_system.control.pulse_dt)`.
- Line 98: conditional branch on `fid<0.80`.

### Key state/data transformations

- Lines 18: computes `sys.magnet` using `sys.magnet=0`.
- Lines 21: computes `sys.isotopes` using `sys.isotopes={'C4','E'}`.
- Lines 24: computes `inter.modes.frqs` using `inter.modes.frqs={0 []}`.
- Lines 27: computes `inter.modes.dispersive` using `inter.modes.dispersive=cell(2,2)`.
- Lines 28: computes `inter.modes.dispersive{1,2}` using `inter.modes.dispersive{1,2}=656.2e3`.
- Lines 31: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `H` using `H=hamiltonian(assume(spin_system,'cavity'))`.
- Lines 42: computes `Cr` using `Cr=operator(spin_system,'C',1); An=operator(spin_system,'A',1)`.
- Lines 43: computes `Lx` using `Lx=operator(spin_system,'Lx',2); Ly=operator(spin_system,'Ly',2)`.
- Lines 44: computes `ops` using `ops={(Cr+An)/2,1i*(Cr-An)/2,Lx,Ly}`.
- Lines 47: computes `rho_init` using `rho_init=state(spin_system,{'BL1','ZL2'},{1,2})`.
- Lines 48: computes `rho_targ` using `rho_targ=state(spin_system,{'BL3','ZL2'},{1,2})`.
- Lines 53: computes `wave_set` using `wave_set=[wave_basis('sine_waves',2,40) wave_basis('cosine_waves',3,40)]'`.
- Lines 56: computes `control.isotopes` using `control.isotopes={'C4','E'}`.
- Lines 57: computes `control.channels` using `control.channels=[1;1;2;2]`.
- Lines 58: computes `control.drifts` using `control.drifts={{H}}`.

## Implementation structure

- GRAPE preparation of a cavity Fock state through a dispersively
- coupled qubit using smooth band-limited drives. The controls are
- expanded in an orthonormal basis of slow sine and cosine waves,
- and the optimisation runs over the expansion coefficients, so
- that the resulting pulses are hardware-friendly smooth envelopes
- rather than free piecewise-constant switches. Compare with the
- piecewise-constant treatment in cavity_fock_grape_a.m; model and
- parameters follow the smooth pulse bosonic GRAPE example of the
- paraqeet package.
- Calculation time: minutes
- Magnet field
- Truncated cavity mode and a qubit

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `state()`, `wave_basis()`, `optimcon()`, `fmaxnewton()`, `pulse()`, `propagator()`, `cumsum()`, `kfigure()`, `kxlabel()`, `kylabel()`, `ktitle()`.
