# examples/quantum_tech/geometric_phase_gate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/geometric_phase_gate.m`
- Signature: `geometric_phase_gate()`
- Total lines: 139

## Purpose

Geometric phase gate between two trapped-ion qubits driven by a state-dependent optical dipole force, as demonstrated by Leibfried et al. (Nature 422, 412 (2003)) with two 9Be+ ions, simulated for a three-ion chain in the same trap: the global force beam pair detuned by delta from the stretch mode pushes the outer ions in opposite directions and leaves the middle ion untouched, so that the outer ions acquire a state-dependent geometric phase while the middle ion is a spectator. Calculation time: seconds

## Physical / mathematical content

- Ion qubits are spin-1/2 particles, the three axial normal modes of the chain (stretch at 6.1 MHz, centre-of-mass at stretch/sqrt(3), Egyptian at sqrt(29/5) times centre-of-mass, James, Appl. Phys. B 66, 181 (1998)) are bosonic modes `V6`, `V3`, `V3`; the spin-dependent force is a static `inter.modes.longitudinal` coupling in the frame rotating with the force drive, and the mode frequencies are laboratory values brought into that frame by `parameters.mode_offset` of the device context under the `spin-phonon` assumption set.
- The stretch mode force constant `kappa=delta/4` closes the phase-space loop at `T=1/delta` with a pi/2 differential geometric phase between the outer ions; the force constants of the other modes follow from the normal mode vectors and the 1/sqrt(frequency) scaling of the zero-point motion.
- Observables are normalised by the overlap with the unit state so that they are true expectation values; the parity of the outer ions after an analysis pi/2 pulse of variable phase reproduces the measurement of Leibfried et al.

## Numerical / algorithmic content

- Basis set `IK-SBS` with `bas.inter_level=[2 3 2]`: boson-boson, spin-boson, and spin-spin coupling graphs traced separately; pure spin-spin correlations up to order two are kept inside the spin-boson subgraphs, which carries the entangling phase. The complete basis for this system has 186624 states.
- The pulse sequence is a local function passed to `device`: a global pi/2 pulse, a `evolution` trajectory over the state-dependent displacement, and a loop of analysis pulses for the parity scan.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-30: Trap parameters, Leibfried et al. (2003); implemented by `stretch_frq=6.1e6;         % stretch mode frequency, Hz`.
- Lines 32-34: Other normal modes of a three-ion chain, James (1998); implemented by `com_frq=stretch_frq/sqrt(3);`.
- Lines 36-37: Force drive frequency; implemented by `drive_frq=stretch_frq+delta;`.
- Lines 39-40: Stretch mode force constant giving a pi/2 geometric phase at T=1/delta; implemented by `kappa=delta/4;`.
- Lines 42-44: Lamb-Dicke scaling of the force constants of the other modes; implemented by `kappa_com=kappa*sqrt(2/3)*sqrt(stretch_frq/com_frq);`.
- Lines 46-48: Ion qubits and the three axial modes of the chain; implemented by `sys.magnet=0;`.
- Lines 50-51: Laboratory frame mode frequencies; implemented by `inter.modes.frqs={[] [] [] stretch_frq com_frq egypt_frq};`.
- Lines 53-62: Spin-dependent forces, Spinach convention Lz(a+a')/sqrt(2); implemented by `inter.modes.longitudinal=cell(6,6);`.
- Lines 64-68: Basis set; implemented by `bas.formalism='sphten-liouv';`.
- Lines 70-72: Spinach housekeeping; implemented by `spin_system=create(sys,inter);`.
- Lines 74-79: Sequence parameters; implemented by `parameters.rho0=state(spin_system,{'ZL1','ZL1','ZL1','BL1','BL1','BL1'},{1,2,3,4,5,6});`.
- Lines 81-82: Run the gate through the device context in the frame of the force drive; implemented by `answer=device(spin_system,@gate_sequence,parameters,'spin-phonon');`.
- Lines 84-90: Validate the spectator ion and the closure of the motional loop; implemented by `if abs(answer.sigma_x(2,end)-answer.sigma_x(2,1))>1e-2`.
- Lines 92-95: Validate the entangling phase: parity amplitude of the outer ions; implemented by `if abs(max(answer.parity)-1)>2e-2`.
- Lines 97-106: Plot the results; implemented by `time_axis=1e6*linspace(0,parameters.gate_time,parameters.npoints);`.
- Lines 110-111: State-dependent displacement sandwiched between two global pi/2 pulses; implemented by `function answer=gate_sequence(spin_system,parameters,H,R,K)`.
- Lines 113-114: Compose the Liouvillian; implemented by `L=H+1i*R+1i*K;`.
- Lines 116-120: Global rotation generators, observables, and the unit state for normalisation; implemented by `Lx=operator(spin_system,'Lx','E'); Ly=operator(spin_system,'Ly','E');`.
- Lines 122-123: First pi/2 pulse puts every ion along +x; implemented by `rho=step(spin_system,Ly,parameters.rho0,pi/2);`.
- Lines 125-129: State-dependent displacement; implemented by `dt=1/parameters.sweep;`.
- Lines 131-136: Parity of the outer ions after an analysis pi/2 pulse of variable phase; implemented by `answer.phases=linspace(0,2*pi,101); answer.parity=zeros(size(answer.phases));`.
