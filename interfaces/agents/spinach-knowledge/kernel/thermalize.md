# kernel/thermalize.m

- Signature: `R=thermalize(spin_system,R,HLSPS,T,rho_eq,method)`

## Purpose and methods

Modifies a relaxation superoperator `R` to drive a chosen stationary state. The two accepted methods modify `R` differently:

- `method='IME'`: requires a nonempty numeric column `rho_eq` and a Liouville-space formalism. In `sphten-liouv`, `U` selects the first basis coordinate; in `zeeman-liouv`, `U` is the vectorised identity. The update is `R=R-kron(U',R*rho_eq)`. For this IME correction to work, the propagated state’s unit-state population must be exactly 1; Spinach cannot check or enforce this requirement. Otherwise the rank-one correction drives towards a mis-scaled equilibrium rather than `rho_eq`.
- `method='dibari'`: requires a nonempty square `HLSPS` and positive real scalar temperature `T`. With `beta=spin_system.tols.hbar/(spin_system.tols.kbol*T)`, the update is `R=R*propagator(spin_system,HLSPS,1i*beta)`. This branch uses the lab-frame Hamiltonian supplied by the caller.

## Input checks

`R` must be numeric and square. Before either method, `norm(R*unit_state(spin_system),2)>1e-10` is rejected as already thermalised. For `dibari`, `norm(HLSPS*unit_state(spin_system),2)<1e-8` is rejected as apparently a commutation superoperator. The source checks these conditions and returns the modified operator; it does not calculate `rho_eq` in the IME branch.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/thermalize.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=thermalize.m)
