# kernel/states/zftrip.m

- Signature: `rho=zftrip(spin_system,ZFS,pops,Z,B,idx)`

## Purpose

Builds a Spinach state from specified zero-field Cartesian triplet populations, then projects it onto the high-field eigenbasis of the combined zero-field-splitting and Zeeman Hamiltonian.

## Inputs and conventions

- `ZFS`: real symmetric 3-by-3 laboratory-frame zero-field-splitting tensor in Hz. The source recommends [`zfs2mat`](../conventions/transforms/zfs2mat.md) for constructing it from `D`, `E`, and molecular Euler angles.
- `pops`: finite, non-negative three-element vector whose sum is 1 within `1e-10`; order is `[pX pY pZ]` for the zero-field X, Y, and Z eigenstates.
- `Z`: real symmetric 3-by-3 laboratory-frame tensor in Hz/Tesla (the code multiplies its Zeeman term by `2*pi`). This is not a gyromagnetic ratio stated in rad/(s*T).
- `B`: real scalar field in Tesla, applied along laboratory Z.
- `idx`: one-based particle index. The source checks `spin_system.comp.isotopes{idx} == 'E3'`; this is isotope metadata for a triplet electron, not a nuclear spin quantum number.

The X/Y/Z labels follow the organic triplet convention `|Dzz|>|Dxx|>|Dyy|`, with `-1/3 < E/D < 0`. Populations expressed in the transition-metal convention `|Dzz|>|Dyy|>|Dxx|`, with `0 < E/D < 1/3`, must swap X and Y before calling. The source notes that at `E=0`, X and Y are degenerate (all three are degenerate when `D=0` too); at `E/D=-1/3`, X and Z have opposite-sign energies with equal magnitudes, so sorting by absolute energy cannot distinguish them. The affected populations must be equal for the result to be meaningful. The cited convention reference is Poole, Farach, and Jackson, *J. Chem. Phys.* 61, 2220 (1974), DOI [10.1063/1.1682294](https://doi.org/10.1063/1.1682294).

## Construction and output

The function forms the spin-1 ZFS Hamiltonian, orders its eigenvectors by increasing absolute eigenvalue to assign the organic X/Y/Z labels, and forms the zero-field density matrix from `pops`. It constructs the Zeeman Hamiltonian using `Z` and `B`, diagonalises the combined Hamiltonian, and removes coherences in that high-field eigenbasis. It then expands the projected operator over nine spherical-tensor labels using `state` calls on particle `idx`, in the order `T0,0`, `T1,+1`, `T1,0`, `T1,-1`, `T2,+2`, `T2,+1`, `T2,0`, `T2,-1`, `T2,-2`. The tensors are scaled by their squared Frobenius norms. No further normalisation of the returned expansion is performed.

The source documents `rho` as a vector in Liouville space. It is assembled in the basis/formalism handled by those `state` calls; this function's own argument checker does not impose a specific `spin_system.bas.formalism`; the called `state` routines handle their own basis checks.

## Links

- Tensor conversion: [`zfs2mat`](../conventions/transforms/zfs2mat.md) and [`axrh2mat`](../conventions/transforms/axrh2mat.md).
- Source: [`kernel/states/zftrip.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/states/zftrip.m).
- Wiki: [`zftrip.m`](https://spindynamics.org/wiki/index.php?title=zftrip.m).
