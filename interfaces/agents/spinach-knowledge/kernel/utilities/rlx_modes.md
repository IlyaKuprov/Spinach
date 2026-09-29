# kernel/utilities/rlx_modes.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rlx_modes.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rlx_modes.m)

## Purpose

Builds the bosonic mode dissipation superoperator: thermalised GKSL dissipators for amplitude damping and pure dephasing of the bosonic modes declared in `inter.modes`, using the amplitude damping rates and pure dephasing rates ingested by `create.m` and the Bose-Einstein thermal occupation numbers computed from the physical mode frequencies and the system temperature.

## Behaviour

- Syntax: `R=rlx_modes(spin_system)`.
- Runs a consistency check (`grumble`) that errors if `spin_system.inter.modes` is missing, if the basis formalism is not `zeeman-liouv` or `sphten-liouv` (bosonic mode dissipators are only available in Liouville space), or if `spin_system.rlx.temperature` is missing.
- Preallocates the output with `mprealloc(spin_system,1)`.
- Locates bosonic modes as components whose type is `C`, `V`, or `T`.
- For each mode, reads the amplitude damping rate `kappa` from `spin_system.inter.modes.damp(k)` and the pure dephasing rate `gphi` from `spin_system.inter.modes.dephase(k)`; if both are zero, the mode is skipped.
- Thermal occupation `nbar` is computed only when `kappa>0` and `spin_system.rlx.temperature>0`:
  - The physical frequency is `carriers(k)+frqs(k)` when `spin_system.inter.modes.carriers(k)>0` (laboratory carrier included where declared), otherwise `abs(spin_system.inter.modes.frqs(k))`.
  - If the physical frequency is below `2*pi*spin_system.tols.inter_cutoff`, the function errors: the thermal occupation of a zero-frequency mode is undefined and the laboratory frame frequency must be supplied.
  - Otherwise `nbar=1/(exp(beta_factor)-1)` with `beta_factor=spin_system.tols.hbar*phys_frq/(spin_system.tols.kbol*spin_system.rlx.temperature)` (Bose-Einstein statistics at the physical frequency), and the used frequency and `nbar` are reported to the user.
  - Without damping or at zero temperature, `nbar=0`. The Spinach convention that zero temperature means the high-temperature limit does not apply to bosonic modes: zero temperature here produces zero thermal occupation numbers.
- Reports `kappa` (s^-1), `gamma_phi` (s^-1) and `nbar` for each dissipative mode.
- Builds left/right ladder superoperators with `operator(spin_system,{'A'|'C'|'N'|'AC'},{k},'left'/'right')`, truncated to the level count of each mode, and accumulates:
  - Cooling dissipator: `kappa*(1+nbar)*(a_left*c_right-0.5*(n_left+n_right))`, i.e. `kappa*(1+nbar)*D[a]` where `D[x]` is the GKSL dissipator of operator `x`.
  - Heating dissipator, added only when `nbar>0`: `kappa*nbar*(c_left*a_right-0.5*(ac_left+ac_right))`, i.e. `kappa*nbar*D[c]`.
  - Pure dephasing dissipator, added only when `gphi>0`: `2*gphi*(n_left*n_right-0.5*(n_left*n_left+n_right*n_right))`, i.e. `2*gamma_phi*D[n]`.

## Inputs and outputs

**Inputs**

- `spin_system` — Spinach spin system description object with bosonic mode information present.

**Outputs**

- `R` — bosonic mode dissipation superoperator.

## References

- Spinach Wiki: [rlx_modes.m](https://spindynamics.org/wiki/index.php?title=rlx_modes.m)
