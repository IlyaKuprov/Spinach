# examples/quantum_tech/optomechanics/optomech_sideband.m

- Signature: `optomech_sideband()`
- Source: [`examples/quantum_tech/optomechanics/optomech_sideband.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/optomechanics/optomech_sideband.m)

## Purpose and model

This is a driven cavity–mechanical-oscillator sideband-transfer simulation, not a spin-defect or EPR calculation. The source declares `sys.isotopes={'C5','V11'}`: a five-level cavity mode and an eleven-level vibrational (phonon) mode. The initial state is cavity level `BL1` and mechanical Fock level `BL3`, i.e. an empty cavity and two phonons in the source's level convention.

In the red-detuned cavity rotating frame, both mode frequencies are set to `10/(2*pi)`; the longitudinal radiation-pressure interaction is `-sqrt(2)/(2*pi)`, coupling cavity occupation to the mechanical coordinate. A coherent cavity drive is added as `2*(C+A)`. The source comments identify this parameter set with the propagation test set of QuantumPropagators.jl. The listed values are dimensionless source units; the frequency and coupling inputs are divided by `2*pi` for the Spinach Hz convention. The model uses `zeeman-hilb` with no basis approximation.

## Propagation and plot

The script constructs the Hamiltonian from these interactions plus the drive, then propagates the density operator with a one-step propagator. It records cavity and mechanical number-operator expectations, `Nc` and `Nm`, on a grid with `dt=0.2`, `250` steps, and dimensionless plotted time from 0 to 50. The plot is the two mode occupations versus that dimensionless time; it is intended to show the sideband-mediated transfer of the two-phonon excitation into the cavity field.

The source includes checks for unit trace, initial mechanical occupation two, and a maximum cavity occupation of at least 1.5. These are source-coded assertions, not stored output values. The script specifies neither damping nor measured device performance.
