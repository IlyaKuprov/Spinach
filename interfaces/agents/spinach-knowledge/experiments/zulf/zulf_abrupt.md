# experiments/zulf/zulf_abrupt.m

- Signature: `fid=zulf_abrupt(spin_system,parameters,H,R,K)`

## Purpose

Models zero-field magnetometry by starting from isotropic thermal equilibrium, propagating through an exponential external-field drop, and acquiring an FID at zero field. The field-drop function is `expdrop`; the subsequent acquisition is simulated with the coupling Hamiltonian and supplied relaxation and kinetics terms.

## Inputs and units

- `spin_system` supplies the spin system, initial magnet field (`spin_system.inter.magnet`), and magnetogyric ratios.
- `parameters.drop_field` is the target field in tesla, non-negative.
- `parameters.drop_time` is the drop duration in seconds and must be positive.
- `parameters.drop_npoints` is a positive integer number of drop discretisation points.
- `parameters.drop_rate` is a positive real scalar in hertz, passed to `expdrop`.
- `parameters.sweep` is the positive acquisition spectral width in hertz; acquisition uses timestep `1/parameters.sweep`.
- `parameters.npoints` is a positive integer acquisition point count.
- `parameters.detection` must be `'uniaxial'` or `'quadrature'`.
- `parameters.flip_angle` is a real numeric scalar in radians, specified for protons; gamma weighting scales the pulse for other nuclei.
- `H`, `R`, and `K` must be numeric matrices of equal dimensions. `R` and `K` enter both field-drop propagation and acquisition. Although accepted and shape-checked, `H` is not used to construct the propagated Hamiltonian: the routine constructs its own field and coupling Hamiltonians from `spin_system`.

## Propagation and readout

The Zeeman and coupling Hamiltonians are generated under the lab-frame `zeeman` and `couplings` assumptions. The Zeeman term is divided by the initial magnet field to normalise it to one tesla. `expdrop` receives the initial field, target field, drop time, drop point count, and drop rate; the propagation timestep is `drop_time/drop_npoints`. Starting from `equilibrium(spin_system)`, each generated field value multiplies the normalised Zeeman term, while the coupling Hamiltonian and `1i*R + 1i*K` are also included in that step.

At the end of the drop, the routine builds gamma-weighted `L+` detection and `Ly` pulse operators using `spin_system.inter.gammas/spin('1H')`, applies the pulse, and acquires using `H_c + 1i*R + 1i*K` with timestep `1/parameters.sweep` and interval count `parameters.npoints-1`. The requested detection mode leaves quadrature data unchanged or returns only the real part for uniaxial detection. The result is the FID on the internally generated gamma-weighted detection state.

## Guardrails and scope

The source checks that `H`, `R`, and `K` are numeric matrices with matching dimensions; sweep and drop rate are positive real scalars; `npoints` and `drop_npoints` are positive integers; target field is a non-negative real scalar; drop time is positive; detection mode is one of the two listed strings; and flip angle is a real numeric scalar. Its header explicitly notes that the supplied offset/Hamiltonian argument is ignored in favor of its internally generated Hamiltonian. Its real-scalar comparisons do not explicitly reject non-finite values.

## Links

- MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/zulf/zulf_abrupt.m
- [Spinach Wiki: zulf_abrupt.m](https://spindynamics.org/wiki/index.php?title=zulf_abrupt.m)
