# kernel/multiprop.m

- Signature: `rho=multiprop(spin_system,P,rho,N)`

## Purpose

Apply the propagator `P` `N` times to `rho` using binary exponentiation, without explicitly constructing `P^N`. The operation depends on `spin_system.bas.formalism`.

## Physical / mathematical content

For `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef`, each active binary power updates the state as `rho=P*rho`. In `zeeman-hilb`, it updates the density matrix as `rho=P*rho*P'`. Thus `N=0` leaves `rho` unchanged, after input validation.

## Numerical / algorithmic content

The routine processes the binary representation of `N`, applying the current propagator power only for set bits and squaring `P` between powers. It cleans the squared propagator using `spin_system.tols.prop_chop` only when higher powers remain to be processed.

Inputs are checked before the zero-step return. The formalism must be one of `zeeman-hilb`, `zeeman-liouv`, `sphten-liouv`, or `zeeman-wavef`; `prop_chop` must be a non-negative real numeric scalar; `P` must be a finite numeric square matrix; and `rho` must be a finite numeric matrix with a row count matching `P`. For `zeeman-hilb`, `rho` must also be square. `N` must be a non-negative real numeric integer scalar; values not represented as integer classes must be finite and no greater than `flintmax`.

## Parameters / inputs

- `spin_system` - Spinach spin-system structure containing `bas.formalism` and `tols.prop_chop`.
- `P` - finite numeric square propagator matrix.
- `rho` - finite numeric state matrix; use state vectors/stacks for the Liouville and wavefunction formalisms, or a density matrix for Hilbert space. Its row count must match `P`; it must be square for `zeeman-hilb`.
- `N` - non-negative integer number of propagator applications, supplied as a real numeric scalar and subject to the source integer-range check.

## Outputs

- `rho` - propagated state; for `N=0`, the validated input is returned unchanged.

## Implementation structure

The implementation validates the formalism, tolerance, matrices, dimensions, and exponent, then performs binary exponentiation. Liouville and wavefunction states are left-multiplied by the active power; Hilbert-space density matrices are transformed on both sides. Propagator powers are cleaned after squaring only when later binary powers are still required.
