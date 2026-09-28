# kernel/rotframe.m

- Signature: `Hr=rotframe(spin_system,H0,H,isotope,order)`

## Purpose

Transforms the laboratory-frame Hamiltonian `H=H0+H1` into a rotating frame referenced to spins of the selected isotope, to the requested perturbation order. The formalism is described in https://doi.org/10.1063/1.4928978.

## Physical / mathematical content

The carrier Hamiltonian `H0` defines the frame. The rotation period depends on the selected isotope's gyromagnetic ratio and the magnetic field; the Hilbert-space period is twice the Liouville-space period.

## Numerical / algorithmic content

The function checks its inputs and calls `intrep` with the period and perturbation order. Its auxiliary-matrix method is faster than the commutator-series and diagonalisation alternatives. Numerical frames are refused for all spins under `nmr` and `cavity`, electrons under `esr`, `deer`, `deer-zz`, and `spin-phonon`, and spin-half nuclei under `qnmr`. Nuclei under electron-only rotating sets and higher-spin nuclei under `qnmr` remain in the laboratory frame; `labframe` retains all spin carriers. The numerical transformation is not implemented for `se_dnp_h+`, `se_dnp_h-`, or `se_dnp_h0`, whose assumptions omit the Zeeman interactions needed to form the laboratory-frame `H0+H1`.

## Parameters / inputs

- `spin_system` — spin system with assumptions set by `assume`; the selected isotope must remain in the laboratory frame under those assumptions.
- `H0` — carrier Hamiltonian defining the rotating frame.
- `H` — laboratory-frame Hamiltonian `H0+H1` to transform.
- `isotope` — character string, such as `'1H'`, identifying the spins used to compute the transformation.
- `order` — perturbation-theory order; may be `inf`.

## Outputs

- `Hr` — rotating-frame Hamiltonian.

## Header notes

Both `H` and `H0` must be Hermitian, and assumption metadata must be present. The function refuses unsupported assumptions before performing the transformation.
