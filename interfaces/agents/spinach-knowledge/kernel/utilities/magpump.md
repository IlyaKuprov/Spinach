# kernel/utilities/magpump.m

## Purpose

Adds phenomenological pumping terms to the relaxation superoperator to enable approximate simulation of CIDNP, PHIP and DNP type effects.

## Behaviour

The function adds pumping as a coupling to the unit state: each substance block of `rate*rho` is added to that block's own unit column `bas.offsets(n)+1` of the relaxation superoperator `R`. The call is

```
R=magpump(spin_system,R,rho,rate)
```

For the pumping to work correctly, `rho` must be the unweighted per-molecule polarisation shape from `coil_state`; the unit populations of the propagated state provide the instantaneous concentrations. A target state may contain components in several substances; each is driven only by its own unit population.

The function is only available in the `sphten-liouv` formalism, and may be called repeatedly if multiple states are pumped.

Consistency checks are enforced by an internal `grumble` function:

- `R` must be a numeric matrix.
- `rho` must be a numeric column vector.
- `rate` must be a finite real scalar.
- `spin_system.bas.formalism` must be `sphten-liouv`.
- Every substance unit coordinate of `rho` must be zero; otherwise an error is raised stating that the unit state cannot be pumped.

## Inputs and outputs

**Inputs**

- `spin_system` — spin system object.
- `R` — relaxation superoperator, from `relaxation()`.
- `rho` — the state to be pumped, from `coil_state()`.
- `rate` — pumping rate, Hz.

**Outputs**

- `R` — modified relaxation superoperator.

## References

- Source: [magpump.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/magpump.m)
- Wiki: [magpump.m](https://spindynamics.org/wiki/index.php?title=magpump.m)
