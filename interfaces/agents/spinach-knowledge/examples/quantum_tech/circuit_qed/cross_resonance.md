# examples/quantum_tech/circuit_qed/cross_resonance.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/circuit_qed/cross_resonance.m`
- Signature: `cross_resonance()`
- Total lines: 143

## Purpose

Cross-resonance gate mechanism between two fixed-frequency transmons in the laboratory frame. The control transmon is driven at the frequency of the target transmon; through the static quadrature coupling, the target then rotates about an axis in its equatorial plane by an angle conditioned on the state of the control. The difference between the two condi- tional rotations is the ZX interaction behind the native two-

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Magnet field; implemented by `sys.magnet=0`.
- Lines 32-33: Two transmons with three levels each; implemented by `sys.isotopes={'T3','T3'}`.
- Lines 35-36: Transmon frequencies and anharmonicities; implemented by `inter.modes.frqs={5.5e9 6.0e9}`.
- Lines 39-40: Quadrature coupling between the transmons; implemented by `inter.modes.exchange=cell(2,2)`.
- Lines 43-44: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 47-48: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 51-52: Laboratory frame drift Hamiltonian; implemented by `H0=full(hamiltonian(assume(spin_system,'labframe')))`.
- Lines 54-55: Quadrature drive operators of the two transmons; implemented by `D1=full(operator(spin_system,'C',1)+operator(spin_system,'A',1))`.
- Lines 58-59: Drive tone parameters, both at the calibrated local oscillator frequency; implemented by `w_lo=2*pi*6.0002e9`.
- Lines 63-64: Flat-top Gaussian envelope parameters; implemented by `t_up=30e-9; t_down=120e-9; t_ramp=15e-9`.
- Lines 66-67: Pauli observables of the target transmon qubit subspace; implemented by `Sx=kron(eye(3),[0 1 0; 1 0 0; 0 0 0])`.
- Lines 71-72: Control transmon in the ground and first excited state; implemented by `unit_mat=eye(9); psi_zero=unit_mat(:,1); psi_one=unit_mat(:,4)`.
- Lines 74-75: Excitation number of the target transmon in the product basis; implemented by `n_targ=kron(ones(3,1),(0:2)')`.
- Lines 77-78: Internal time stepping of the source calculation; implemented by `nsub=15000; dts=1e-11; nrec=75`.
- Lines 80-81: Preallocate conditional Bloch trajectories; implemented by `bloch=zeros(3,nsub/nrec+1,2)`.
- Lines 84-85: Propagate on the internal midpoint grid; implemented by `for k=1:nsub`.
- Lines 87-88: Flat-top Gaussian envelope at the midpoint; implemented by `tmid=(k-0.5)*dts`.
- Lines 91-92: IQ mixer signals of the two tones; implemented by `sig_one=amp_one*env*cos(w_lo*tmid-phi_one)`.

### Control flow inferred from the code

- Line 85: `for` loop over `k=1:nsub`.
- Line 100: conditional branch on `mod(k,nrec)==0`.
- Line 114: conditional branch on `(abs(norm(psi_zero)-1)>1e-6)||(abs(norm(psi_one)-1)>1e-6)`.
- Line 119: conditional branch on `~((bloch(2,end,1)>0.5)&&(bloch(3,end,1)>0.2))`.
- Line 124: conditional branch on `~((bloch(2,end,2)>0)&&(bloch(3,end,2)<-0.5))`.

### Key state/data transformations

- Lines 30: computes `sys.magnet` using `sys.magnet=0`.
- Lines 33: computes `sys.isotopes` using `sys.isotopes={'T3','T3'}`.
- Lines 36: computes `inter.modes.frqs` using `inter.modes.frqs={5.5e9 6.0e9}`.
- Lines 37: computes `inter.modes.anharms` using `inter.modes.anharms={-240e6 -200e6}`.
- Lines 40: computes `inter.modes.exchange` using `inter.modes.exchange=cell(2,2)`.
- Lines 41: computes `inter.modes.exchange{1,2}` using `inter.modes.exchange{1,2}=25e6`.
- Lines 44: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 45: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 48: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 52: computes `H0` using `H0=full(hamiltonian(assume(spin_system,'labframe')))`.
- Lines 55: computes `D1` using `D1=full(operator(spin_system,'C',1)+operator(spin_system,'A',1))`.
- Lines 56: computes `D2` using `D2=full(operator(spin_system,'C',2)+operator(spin_system,'A',2))`.
- Lines 59: computes `w_lo` using `w_lo=2*pi*6.0002e9`.
- Lines 60: computes `amp_one` using `amp_one=5.969026041820607e8; phi_one=0`.
- Lines 61: computes `amp_two` using `amp_two=2.883982055995430e7; phi_two=-0.39720756`.
- Lines 64: computes `t_up` using `t_up=30e-9; t_down=120e-9; t_ramp=15e-9`.
- Lines 67: computes `Sx` using `Sx=kron(eye(3),[0 1 0; 1 0 0; 0 0 0])`.
- Lines 68: computes `Sy` using `Sy=kron(eye(3),[0 -1i 0; 1i 0 0; 0 0 0])`.

## Implementation structure

- Cross-resonance gate mechanism between two fixed-frequency
- transmons in the laboratory frame. The control transmon is
- driven at the frequency of the target transmon; through the
- static quadrature coupling, the target then rotates about an
- axis in its equatorial plane by an angle conditioned on the
- state of the control. The difference between the two condi-
- tional rotations is the ZX interaction behind the native two-
- qubit gate of fixed-frequency superconducting processors; at
- these parameters the unconditional part is the larger of the
- two, and an echo sequence would be needed to remove it. A weak
- direct tone on the target sets the gate phases. The simulation
- runs in the laboratory frame, but the recorded kets are rotated

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `unit_mat()`, `bloch()`, `erf()`, `kfigure()`, `scale_figure()`, `subplot()`, `squeeze()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
