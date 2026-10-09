# kernel/homospoil.m

- Signature: `rho=homospoil(spin_system,rho,zqc_flag)`

## Purpose

Projects the supplied state onto the components retained by the selected homospoil approximation. It filters state coefficients; it does not change a Hamiltonian or apply a phase evolution. Chemical-shift offsets are not part of the frequency test.

## Inputs and shape

- `rho` must be numeric; the source documentation describes a state vector or horizontal stack of states. In the Liouville paths, the implementation folds the spin-space basis rows against all remaining columns for filtering and restores the original input shape.
- `zqc_flag` must be the character value `'keep'` or `'destroy'`. It affects the `sphten-liouv` path; the two Zeeman paths ignore it.

## Formalisms and retained components

- In `sphten-liouv`, each substance’s local basis labels are converted to projection indices `M`, and its base frequencies are selected through `chem.parts{n}`. The mask acts only within that block. With `'keep'`, a row survives when `abs(sum(basefrqs .* M,2)) <= 1e-6`; the signed carrier-frequency-weighted sum is used, so contributions from different spins can cancel. With `'destroy'`, only rows with `sum(abs(M),2) == 0` survive, i.e. longitudinal components with zero coherence order on every spin. The frequency test uses `spin_system.inter.basefrqs` directly, with no conversion in this function; the `1e-6` tolerance is in the units of those entries.
- In `zeeman-hilb`, the implementation takes `diag(rho)` and returns a matrix with that diagonal and zero off-diagonal elements, for either flag.
- In `zeeman-liouv`, the spin-space diagonal of every folded Liouville block is retained for either flag; off-diagonal spin-space elements are zeroed and the original state shape is restored.
- Fokker–Planck direct-product dimensions are supported in the Liouville-space formalisms.

The implementation reports a warning if the retained state has 1-norm below `1e-10`.

## Source links

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/homospoil.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=homospoil.m)

Multi-substance Zeeman filtering is unsupported and raises `Spinach:homospoil:segmentedZeeman` before constructing a tensor-product channel or applying a diagonal mask. Single-substance Zeeman and segmented spherical-tensor paths remain available.
