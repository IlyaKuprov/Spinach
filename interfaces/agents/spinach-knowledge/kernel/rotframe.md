# kernel/rotframe.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/rotframe.m`
- Signature: `Hr=rotframe(spin_system,H0,H,isotope,order)`
- Total lines: 115

## Purpose

Rotating frame transformation with respect to specified spins to specified order in perturbation theory, using the formalism described in https://doi.org/10.1063/1.4928978 Syntax: Hr=rotframe(spin_system,H0,H,isotope,order)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

Numerical frames reject all spins under `nmr` and `cavity`; electrons under `esr`, `deer`, `deer-zz`, and `spin-phonon`; and spin-half nuclei under `qnmr`. Nuclei under electron-only rotating sets and higher-spin nuclei under `qnmr` remain in the laboratory frame and may be transformed. `labframe` retains all spin carriers. Numerical frames are not implemented for `se_dnp_h+`, `se_dnp_h-`, or `se_dnp_h0`: `assume` omits all Zeeman interactions from these solid-effect components, so they do not supply the laboratory Hamiltonian H0+H1 required by this transformation. This refusal does not alter the component construction or the solid-effect experiment.

## Parameters / inputs

- spin_system - spin system with assumptions set by `assume`; the selected isotope must still be in the laboratory frame.
- H0 -carrier Hamiltonian with respect to which the
- rotating frame transformation is to be done
- H -laboratory frame Hamiltonian H0+H1 that is to
- be transformed into the rotating frame
- isotope -string, such as '1H', specifying the spins
- with respect to which the transformation is
- being computed
- order -perturbation theory order in the rotating
- frame transformation, this may be inf

## Outputs

- Hr -rotating frame Hamiltonian
- Notes: the auxiliary matrix method is massively faster than
- either commutator series or diagonalisation.

## Header notes

The auxiliary-matrix method used by `intrep` is faster than the commutator-series and diagonalisation alternatives described in the source. The rotation period follows the selected isotope's gyromagnetic ratio and field, with the Hilbert-space period twice the Liouville-space period. `H` and `H0` must be Hermitian; missing assumption metadata is refused before the numerical transformation.
