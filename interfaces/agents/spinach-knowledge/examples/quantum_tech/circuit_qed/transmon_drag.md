# examples/quantum_tech/circuit_qed/transmon_drag.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/circuit_qed/transmon_drag.m`
- Signature: `transmon_drag()`
- Total lines: 129

## Purpose

DRAG correction of a resonant Gaussian pulse on a three-level Duffing transmon in the laboratory frame. The derivative of the envelope, scaled by the DRAG detuning parameter, is fed into the quadrature channel of the IQ mixer; this suppresses the leakage into the second excited state and improves the fidelity of the 90 degree rotation on the qubit subspace. A further improvement comes from numerical optimisation of t

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Magnet field; implemented by `sys.magnet=0`.
- Lines 29-30: Particle specification; implemented by `sys.isotopes={'T3'}`.
- Lines 32-33: Transmon frequency and anharmonicity; implemented by `frq_01=4.8e9; anharm=-200e6`.
- Lines 37-38: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 41-42: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 45-46: Laboratory frame drift Hamiltonian; implemented by `H0=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 48-49: Transmon quadrature drive operator; implemented by `Cr=operator(spin_system,'C',1); An=operator(spin_system,'A',1)`.
- Lines 52-53: Pulse duration and Gaussian envelope width; implemented by `t_dur=20e-9; sigma=t_dur/8`.
- Lines 55-56: Analytical DRAG detuning parameter for a Duffing transmon; implemented by `drag_det=2*(2*pi)*anharm`.
- Lines 58-59: Rows: plain Gaussian, analytical DRAG, optimised DRAG; implemented by `pulse_params=[3.000000000000000e8 drag_det 2*pi*frq_01 0 0`.
- Lines 63-64: Target gate is a 90 degree rotation on the qubit subspace; implemented by `target=[cos(pi/4) -1i*sin(pi/4) 0; -1i*sin(pi/4) cos(pi/4) 0; 0 0 1]`.
- Lines 66-67: Internal time stepping of the source calculation; implemented by `nsub=2000; dts=1e-11`.
- Lines 69-70: Preallocate fidelities and population trajectories; implemented by `fids=zeros(1,3); pops=zeros(3,nsub+1,3)`.
- Lines 72-73: Loop over the three pulses; implemented by `for m=1:3`.
- Lines 75-76: Extract mixer parameters; implemented by `amp=pulse_params(m,1); delta=pulse_params(m,2)`.
- Lines 80-81: Propagate on the internal midpoint grid; implemented by `U=eye(3); psi=[1; 0; 0]; pops(:,1,m)=abs(psi).^2`.
- Lines 84-85: Gaussian envelope and its derivative at the midpoint; implemented by `tmid=(k-0.5)*dts`.
- Lines 89-90: IQ mixer signal with the DRAG quadrature; implemented by `sig=env*cos(w_lo*tmid-phi)-drag*(denv/delta)*sin(w_lo*tmid-phi)`.

### Control flow inferred from the code

- Line 73: `for` loop over `m=1:3`.
- Line 82: `for` loop over `k=1:nsub`.
- Line 107: conditional branch on `max(abs(fids-[0.969588 0.973797 0.994239]))>1e-3`.
- Line 112: conditional branch on `~(pops(3,end,2)<pops(3,end,1)/2)`.

### Key state/data transformations

- Lines 27: computes `sys.magnet` using `sys.magnet=0`.
- Lines 30: computes `sys.isotopes` using `sys.isotopes={'T3'}`.
- Lines 33: computes `frq_01` using `frq_01=4.8e9; anharm=-200e6`.
- Lines 34: computes `inter.modes.frqs` using `inter.modes.frqs={frq_01}`.
- Lines 35: computes `inter.modes.anharms` using `inter.modes.anharms={anharm}`.
- Lines 38: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 39: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 42: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 46: computes `H0` using `H0=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 49: computes `Cr` using `Cr=operator(spin_system,'C',1); An=operator(spin_system,'A',1)`.
- Lines 50: computes `Dx` using `Dx=full(Cr+An); H0=full(H0)`.
- Lines 53: computes `t_dur` using `t_dur=20e-9; sigma=t_dur/8`.
- Lines 56: computes `drag_det` using `drag_det=2*(2*pi)*anharm`.
- Lines 59: computes `pulse_params` using `pulse_params=[3.000000000000000e8 drag_det 2*pi*frq_01 0 0`.
- Lines 64: computes `target` using `target=[cos(pi/4) -1i*sin(pi/4) 0; -1i*sin(pi/4) cos(pi/4) 0; 0 0 1]`.
- Lines 67: computes `nsub` using `nsub=2000; dts=1e-11`.
- Lines 70: computes `fids` using `fids=zeros(1,3); pops=zeros(3,nsub+1,3)`.
- Lines 76: computes `amp` using `amp=pulse_params(m,1); delta=pulse_params(m,2)`.

## Implementation structure

- DRAG correction of a resonant Gaussian pulse on a three-level
- Duffing transmon in the laboratory frame. The derivative of the
- envelope, scaled by the DRAG detuning parameter, is fed into the
- quadrature channel of the IQ mixer; this suppresses the leakage
- into the second excited state and improves the fidelity of the
- 90 degree rotation on the qubit subspace. A further improvement
- comes from numerical optimisation of the mixer parameters. The
- rows of pulse_params are the plain Gaussian, the analytically
- corrected DRAG, and the numerically optimised DRAG pulse; the
- columns are the envelope amplitude, the DRAG detuning, the lo-
- cal oscillator frequency, the mixer phase, and the DRAG quad-
- rature switch. The propagator is taken into the frame of the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `pulse_params()`, `pops()`, `fids()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
