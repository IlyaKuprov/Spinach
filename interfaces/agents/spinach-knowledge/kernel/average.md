# kernel/average.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/average.m`
- Signature: `H=average(spin_system,Hp,H0,Hm,omega,theory)`
- Total lines: 209

## Purpose

Average Hamiltonian theories under Zeeman interaction rotating frame transformations. Syntax: H=average(spin_system,Hp,H0,Hm,omega,theory)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 57-58: Check consistency; implemented by `grumble(Hp,H0,Hm,omega,theory)`.
- Lines 60-61: Print diagnostics; implemented by `report(spin_system,['H+ 1-norm: ' num2str(norm(Hp,1))])`.
- Lines 66-67: Run the averaging; implemented by `switch theory`.
- Lines 150-151: Discretise the period of the rotating frame, 16 intervals are enough; implemented by `nslices=16; time_grid=linspace(0,2*pi/omega,nslices+1); time_step=time_grid(2)`.
- Lines 153-154: Compute the product integral using 4th order Lie quadrature; implemented by `P=eye(size(H0)); report(spin_system,'computing product integral ')`.
- Lines 163-164: Compute matrix logarithm (slow, needs work); implemented by `report(spin_system,'computing matrix logarithm ')`.
- Lines 169-170: Complain and bomb out; implemented by `error('unknown theory specification.')`.
- Lines 174-175: Clean up the result and force the matrix into sparse format; implemented by `H=clean_up(spin_system,H,spin_system.tols.liouv_zero); H=sparse(H)`.
- Lines 177-180: Report to the user; implemented by `report(spin_system,['average Hamiltonian dimension ' num2str(size(H,1)) ', nnz ' num2str(nnz(H)) ', density ' num2str(100*nnz(H)/numel(H)) '%, 1-norm ' num2str(norm(H,1)…`.

### Control flow inferred from the code

- Line 67: dispatches on `theory`; cases `'kb_first_order'`, `'kb_second_order'`, `'kb_third_order'`, `'ah_first_order'`, `'ah_second_order'`, `'ah_third_order'`, `'matrix_log'`.
- Line 155: `parfor` loop over `n=1:nslices`.

### Key state/data transformations

- Lines 71: computes `H1` using `H1=+(1/omega^1)*(Hp*Hm-Hm*Hp); H=H0+H1`.
- Lines 78: computes `H2` using `H2=-(1/omega^2)*(Hp*(Hm*H0-H0*Hm)+Hm*(Hp*H0-H0*Hp))`.
- Lines 79: computes `H` using `H=H0+H1+H2`.
- Lines 88-91: computes `H3` using `H3=+(1/omega^3)*(0.5*(Hm*Hm*Hp*Hp-Hp*Hp*Hm*Hm)+ (Hp*Hm+Hm*Hp)*(Hp*Hm-Hm*Hp)+ Hp*(H0*(Hm*H0-H0*Hm)-(Hm*H0-H0*Hm)*H0)- Hm*(H0*(Hp*H0-H0*Hp)-(Hp*H0-H0*Hp)*H0))`.
- Lines 151: computes `nslices` using `nslices=16; time_grid=linspace(0,2*pi/omega,nslices+1); time_step=time_grid(2)`.
- Lines 154: computes `P` using `P=eye(size(H0)); report(spin_system,'computing product integral ')`.
- Lines 157: computes `HL` using `HL=H0+exp(+1i*omega*(time_grid(n)+0.0*time_step))*Hp+exp(-1i*omega*(time_grid(n)+0.0*time_step))*Hm`.
- Lines 158: computes `HM` using `HM=H0+exp(+1i*omega*(time_grid(n)+0.5*time_step))*Hp+exp(-1i*omega*(time_grid(n)+0.5*time_step))*Hm`.
- Lines 159: computes `HR` using `HR=H0+exp(+1i*omega*(time_grid(n)+1.0*time_step))*Hp+exp(-1i*omega*(time_grid(n)+1.0*time_step))*Hm`.

### Local helper functions

- Line 185: `grumble()` — `function grumble(Hp,H0,Hm,omega,theory)`.
  - Representative operation: `if (~isnumeric(Hp))||(~isnumeric(H0))||(~isnumeric(Hm))|| (~ismatrix(Hp))||(~ismatrix(H0))||(~ismatrix(Hm))`.
  - Representative operation: `(~ismatrix(Hp))||(~ismatrix(H0))||(~ismatrix(Hm))`.

## Parameters / inputs

- Hp -the part of the rotating frame Hamiltonian that has positive
- frequency +omega under the rotating frame transformation
- H0 -the part of the rotating frame Hamiltonian that has zero
- frequency under the rotating frame transformation
- Hm -the part of the rotating frame Hamiltonian that has negative
- frequency -omega under the rotating frame transformation
- omega -the frequency of the rotating frame transformation, rad/s
- theory -the level of the average Hamiltonian theory:
- 'ah_first_order' -first order in Waugh theory
- 'ah_second_order' -second order in Waugh theory
- 'ah_third_order' -third order in Waugh theory
- 'matrix_log' -exact algorithm (very expensive,
- uses dense matrix algebra)
- 'kb_first_order' -first order in Krylov-Bogolyubov
- theory (DNP experiments only)
- 'kb_second_order' -second order in Krylov-Bogolyubov
- theory (DNP experiments only)
- 'kb_third_order' -third order in Krylov-Bogolyubov
- theory (DNP experiments only)

## Outputs

- H -average Hamiltonian
- Note: Krylov-Bogolyubov averging theory as applied to DNP systems is
- described in detail here:

## Implementation structure

- Average Hamiltonian theories under Zeeman interaction rotating frame
- transformations. Syntax:
- H=average(spin_system,Hp,H0,Hm,omega,theory)
- Hp - the part of the rotating frame Hamiltonian that has positive
- frequency +omega under the rotating frame transformation
- H0 - the part of the rotating frame Hamiltonian that has zero
- frequency under the rotating frame transformation
- Hm - the part of the rotating frame Hamiltonian that has negative
- frequency -omega under the rotating frame transformation
- omega - the frequency of the rotating frame transformation, rad/s
- theory - the level of the average Hamiltonian theory:
- 'ah_first_order' -first order in Waugh theory

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `num2str()`, `time_grid()`, `propagator()`, `isergen()`, `logm()`, `clean_up()`, `nnz()`, `issparse()`, `ismatrix()`, `any()`, `isscalar()`, `ischar()`.
