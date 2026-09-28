# kernel/correlation.m

- Signature: `rho=correlation(spin_system,rho,orders,spins)`

## Purpose

Keeps only the requested correlation orders in a state vector or a stack of state vectors. In Zeeman-Hilbert formalism, the input and output may instead be density matrices or stacks of density matrices. The selection can serve as an analytical alternative to complicated phase cycles.

## Physical / mathematical content

Correlation order is evaluated over the selected spins. In the spherical-tensor Liouville formalism, the function counts selected spins with nonzero basis components and retains basis states whose count is in `orders`. In Zeeman Liouville and Zeeman-Hilbert formalisms, it projects onto the requested orders using per-spin identity-component channels; the function documentation notes that correlation order is not diagonal in the Zeeman basis.

## Numerical / algorithmic content

- In `sphten-liouv`, a mask of basis states with the requested correlation orders is applied to every column of the input.
- In `zeeman-liouv` and `zeeman-hilb`, per-spin identity-channel projections are combined using roots-of-unity samples and discrete Fourier weights to select the requested orders. Zeeman-Hilbert density matrices are expanded into Liouville space for filtering and then folded back.
- The input dimensions are restored after filtering. If the resulting one-norm is below `1e-10`, the function reports a warning that magnetization appears to have been destroyed.

The supported formalisms are `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb`. Fokker–Planck direct products are supported in the Liouville-space formalisms.

## Parameters / inputs

- `spin_system` — Spinach system structure with basis and formalism information.
- `rho` — numeric state vector or horizontal stack; in `zeeman-hilb`, a density matrix or horizontal stack of density matrices.
- `orders` — vector of non-negative integer correlation orders to retain (called `correlation_orders` in the parameter description).
- `spins` — optional spin selection: `'all'` by default, an isotope label such as `'1H'` or `'13C'`, or a vector of spin numbers.

## Outputs

- `rho` — filtered state vector(s) or density matrix/matrices, with unrequested correlation orders zeroed; the input dimensions are preserved.

## Reference

- [Spin Dynamics Wiki: correlation.m](https://spindynamics.org/wiki/index.php?title=correlation.m)
