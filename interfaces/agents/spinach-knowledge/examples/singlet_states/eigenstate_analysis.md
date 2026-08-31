# examples/singlet_states/eigenstate_analysis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/eigenstate_analysis.m`
- Signature: `eigenstate_analysis()`
- Total lines: 100

## Purpose

Stationary state analysis for the spin system of allyl pyruvate, finding out which component of the singlet state commutes with the drift Hamiltonian.

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Get the spin system from Anu's fits; implemented by `[sys,inter]=allyl_pyruvate({'1H','13C'})`.
- Lines 12-13: Set the magnet; implemented by `sys.magnet=14.1`.
- Lines 15-16: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 18-19: Pick out the required 13C isotopomer; implemented by `subsystems=dilute(spin_system,'13C',1)`.
- Lines 22-23: Generate the basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 27-28: Get isotropic Hamiltonian; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 31-32: Tidy up rounding errors; implemented by `H=(H+ctranspose(H))/2`.
- Lines 34-35: Get the singlet state; implemented by `rho=singlet(spin_system,3,4)`.
- Lines 37-39: Report the norm; implemented by `report(spin_system,['Singlet state norm: ' num2str(norm(full(rho),2))])`.
- Lines 41-42: Remove the part that does not commute with H; implemented by `[EvH,~]=eig(H); rho_inv=remncomm(rho,EvH)`.
- Lines 44-45: Remove the unit part; implemented by `rho_inv=remtrace(rho_inv)`.
- Lines 47-49: Report the norm; implemented by `report(spin_system,[' of which commutes with H: ' num2str(norm(full(rho_inv),2))])`.
- Lines 51-52: Project out the ZZ component; implemented by `zz_state=state(spin_system,{'Lz','Lz'},{3 4})`.
- Lines 56-58: Report the norm; implemented by `report(spin_system,[' of which is not ZZ: ' num2str(norm(full(rho_inv),2))])`.
- Lines 60-63: % Apply offset and spin-lock, repeat process; implemented by `parameters.offset=2850`.
- Lines 62-63: Hamiltonian offset; implemented by `parameters.offset=2850`.
- Lines 67-68: Spin-lock term at 1 kHz; implemented by `H=H+2*pi*1000*operator(spin_system,'Lx','1H')`.

### Key state/data transformations

- Lines 10: computes `[sys,inter]` using `[sys,inter]=allyl_pyruvate({'1H','13C'})`.
- Lines 13: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 16: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 19: computes `subsystems` using `subsystems=dilute(spin_system,'13C',1)`.
- Lines 23: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `H` using `H=hamiltonian(spin_system); H=full(H)`.
- Lines 35: computes `rho` using `rho=singlet(spin_system,3,4)`.
- Lines 42: computes `[EvH,~]` using `[EvH,~]=eig(H); rho_inv=remncomm(rho,EvH)`.
- Lines 45: computes `rho_inv` using `rho_inv=remtrace(rho_inv)`.
- Lines 52: computes `zz_state` using `zz_state=state(spin_system,{'Lz','Lz'},{3 4})`.
- Lines 63: computes `parameters.offset` using `parameters.offset=2850`.
- Lines 64: computes `parameters.spins` using `parameters.spins={'1H'}`.

## Implementation structure

- Stationary state analysis for the spin system of allyl pyruvate,
- finding out which component of the singlet state commutes with
- the drift Hamiltonian.
- Get the spin system from Anu's fits
- Set the magnet
- Spinach housekeeping
- Pick out the required 13C isotopomer
- Generate the basis
- Get isotropic Hamiltonian
- Tidy up rounding errors
- Get the singlet state
- Report the norm

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `allyl_pyruvate()`, `create()`, `dilute()`, `basis()`, `assume()`, `hamiltonian()`, `ctranspose()`, `singlet()`, `report()`, `num2str()`, `remncomm()`, `remtrace()`, `state()`, `frqoffset()`, `operator()`.
