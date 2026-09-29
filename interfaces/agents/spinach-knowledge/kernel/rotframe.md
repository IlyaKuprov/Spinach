# kernel/rotframe.m

- Signature: `Hr=rotframe(spin_system,H0,H,isotope,order)`

## Purpose

Transforms the laboratory-frame Hamiltonian `H=H0+H1` into a rotating frame referenced to the selected isotope, to the requested perturbation order. The source cites the formalism in [doi:10.1063/1.4928978](https://doi.org/10.1063/1.4928978).

## Physical / mathematical content

`H0` is the carrier Hamiltonian defining the frame; `H` is the laboratory-frame Hamiltonian `H0+H1`. The selected isotope is a character string such as `'1H'`, and `order` is the perturbation-theory order, which may be `inf`. The source computes the period from the isotope gyromagnetic ratio and field: `T=-2*pi/(spin(isotope)*spin_system.inter.magnet)` for Liouville-space formalisms and `T=-4*pi/(spin(isotope)*spin_system.inter.magnet)` for Hilbert-space formalisms.

## Numerical / algorithmic content

After validating the inputs, the wrapper passes `H0`, `H`, the selected period, and `order` to `intrep`. Both Hamiltonians must be Hermitian, assumption metadata must be present, and the selected isotope must still be in the laboratory frame under those assumptions.

Numerical frames are refused for all spins under `nmr` and `cavity`, electrons under `esr`, `deer`, `deer-zz`, and `spin-phonon`, and spin-1/2 nuclei under `qnmr`. Nuclei under electron-only rotating sets and higher-spin nuclei under `qnmr` remain in the laboratory frame. Carrier-free solid-effect components `se_dnp_h+`, `se_dnp_h-`, and `se_dnp_h0` are not supported as numerical frames because the source identifies them as not being laboratory Hamiltonians of the form `H0+H1`.

## Parameters / inputs

- `spin_system` — spin system with assumptions set by `assume()`.
- `H0` — carrier Hamiltonian defining the rotating frame.
- `H` — laboratory-frame Hamiltonian to transform.
- `isotope` — character string identifying spins used to compute the transformation, for example `'1H'`.
- `order` — perturbation-theory order; may be `inf`.

## Outputs

- `Hr` — rotating-frame Hamiltonian.

## Source

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/rotframe.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=rotframe.m)
