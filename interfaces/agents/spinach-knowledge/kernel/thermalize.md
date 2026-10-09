# kernel/thermalize.m

- Signature: `R=thermalize(spin_system,R,HLSPS,T,rho_eq,method)`

## Purpose and methods

Modifies a relaxation superoperator `R` to drive a chosen stationary state. The two accepted methods modify `R` differently:

- `method='IME'`: requires a nonempty numeric column `rho_eq` and a Liouville-space formalism. The correction is applied independently within each substance block: `R_n=R_n-(R_n*rho_eq_n)*U_n'`. In `sphten-liouv`, `U_n` selects that block’s first coordinate; in `zeeman-liouv`, it is the vectorised local identity. The supplied target blocks are unweighted unit-concentration shapes, requested from `equilibrium` with a copy whose `chem.concs` entries are all one. Acting on a concentration-weighted state therefore supplies the target concentration through its unit coordinate, not twice. No concentration division or dynamic block renormalisation occurs. A spherical-tensor target whose unit coordinates differ from one, or a Zeeman-Liouville target whose local trace differs from one, is rejected with `Spinach:thermalize:targetConcentration`. This preserves the direct-sum structure rather than connecting different unit coordinates. Nonzero cross-substance blocks of `R` are unsupported and raise `Spinach:thermalize:crossSubstanceRelaxation` before correction, even if they annihilate the unit state.
- `method='dibari'`: requires a nonempty square `HLSPS` and positive real scalar temperature `T`. With `beta=spin_system.tols.hbar/(spin_system.tols.kbol*T)`, the update is `R=R*propagator(spin_system,HLSPS,1i*beta)`. This branch uses the lab-frame Hamiltonian supplied by the caller.

## Hilbert matrix IME

In `zeeman-hilb`, IME accepts a numeric relaxation map on the direct sum of vectorised local matrices, with dimension `sum(D_n^2)`, and a block-diagonal target matrix with unit trace in every block. The same Zeeman-Liouville IME correction is applied, and the returned handle `R(t,rho)` evaluates the block-diagonal matrix derivative through `hilb_action`. It is a matrix RHS, not a Hamiltonian or a `step` generator. This explicit-map interface does not alter the stock Hilbert `relaxation` constructor or its theory restrictions. Invalid map dimensions and non-block-diagonal targets have named rejections; unit target traces and cross-substance relaxation are checked by the shared Liouville path.

## Input checks

`R` must be numeric and square. The guard unit is geometric, constructed with concentrations set to one, so zero-population blocks remain checked. Before either method, `norm(R*unit,2)>1e-10` is rejected as already thermalised. For `dibari`, `norm(HLSPS*unit,2)<1e-8` is rejected as apparently a commutation superoperator. The source checks these conditions and returns the modified operator; it does not calculate `rho_eq` in the IME branch.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/thermalize.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=thermalize.m)

Wavefunction thermalisation is explicitly rejected before constructing any unit state, with `Spinach:thermalize:wavefunction`.
