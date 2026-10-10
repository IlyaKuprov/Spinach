# kernel/correlation.m

Source: [kernel/correlation.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/correlation.m)
Wiki: [Spin Dynamics Wiki: correlation.m](https://spindynamics.org/wiki/index.php?title=correlation.m)

- Signature: `rho=correlation(spin_system,rho,correlation_orders,spins)`

## Purpose

Keeps the requested spin-correlation orders and zeros the other components in the supplied state. The header describes this as an analytical alternative to complicated phase cycles. The implementation uses a basis mask in spherical-tensor Liouville formalism and a discrete Fourier projector in Zeeman Liouville and Hilbert formalisms.

## Inputs and output

- `spin_system` — Spinach system object; this routine reads its basis formalism and basis and the component isotope, multiplicity, and spin-count metadata.
- `rho` — numeric state vector or horizontal stack of state vectors. In `zeeman-hilb`, it may instead be a density matrix or horizontal stack of density matrices. The original shape is restored on output.
- `correlation_orders` — numeric vector of non-negative integer orders to retain. The header describes a row vector; executable validation accepts a numeric vector without requiring row orientation. Orders outside `0:numel(spins selected)` have no corresponding Zeeman Fourier component.
- `spins` — optional selector: `'all'`, an isotope label such as `'1H'` or `'13C'`, or a vector of valid one-based spin indices. If omitted, it defaults to `'all'`.
- Output `rho` has the input dimensions, with components at unrequested orders set to zero.

## How the selection is applied

The code takes `bas.offsets(end)` as the spin-space dimension, squares it for `zeeman-hilb`, and treats the remaining flattened columns as the space dimension. It reshapes `rho` for filtering and restores its original dimensions afterward.

For `sphten-liouv`, the routine counts nonzero entries in each local descriptor’s columns belonging to the selected global spins, mapped through `chem.parts{n}`. The resulting masks are placed at `bas.offsets(n)`. It retains basis rows whose selected-spin correlation order is requested and zeros all other rows across the input columns.

For `zeeman-liouv` and `zeeman-hilb`, it builds sparse identity-component channels for the selected spins, samples the generating operation at roots of unity, and combines the samples with discrete Fourier weights for the requested orders from zero through the number of selected spins. In Hilbert formalism the density matrix is processed through the corresponding squared spin-space dimension; the code restores the original input shape. The header notes that correlation order is not diagonal in the Zeeman basis and describes the Hilbert-space density-matrix handling as a Liouville-space stretch/filter/fold operation; the executable path performs the reshape and projection directly.

The routine accepts only `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb` formalisms. It checks that `rho` is numeric, that orders are a numeric vector of non-negative integers, and that `spins` is either a numeric vector of valid integer indices or a character selector equal to `'all'` or an isotope label present in the system. After filtering it warns when `norm(rho,1)<1e-10`, reporting that all magnetisation appears to have been destroyed.

## Units and examples

Correlation orders and numeric spin selectors are integer indices; no physical units are specified. The source gives the call syntax and selector examples `'1H'`, `'13C'`, and `'all'`, but no worked numerical example. The header also notes support for Fokker–Planck direct-product spaces in Liouville formalisms; the executable switch handles those formalisms through its basis-space dimensions rather than a separate Fokker–Planck branch.

## Signature clarification

The header names the third argument `correlation_orders`; the executable declaration calls it `orders`. This is a parameter-name clarification only.

Multi-substance Zeeman filtering is unsupported and raises `Spinach:correlation:segmentedZeeman` before constructing a tensor-product channel or applying a diagonal mask. Single-substance Zeeman and segmented spherical-tensor paths remain available.

Legacy global `bas.basis` matrices and `bas.irrep` fields are rejected at this entry point with named errors pointing to per-substance `bas.basis{n}`/`bas.offsets` and `bas.sym_fact(n)` symmetry data. Compiled structures remain ordinary MATLAB structs; arbitrary external dot reads are not intercepted.
