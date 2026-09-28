# kernel/utilities/phan2fpl.m

- Signature: `rho=phan2fpl(phan,rho)`

## Purpose

Projects a spatial intensity distribution into Fokker-Planck space as the image painted by the supplied spin state.

## Parameters / inputs

- `phan` — phantom containing the spatial distribution of the amplitude of the specified spin state; it must be numeric, real, and one-, two-, or three-dimensional.
- `rho` — numeric, single-column Liouville-space state vector.

## Output

- `rho` — the corresponding Fokker-Planck-space state vector.

## Implementation structure

After checking the inputs, the function computes `kron(phan(:),rho)`, flattening the phantom into a column before taking the Kronecker product.

## Reference

- <https://spindynamics.org/wiki/index.php?title=phan2fpl.m>
