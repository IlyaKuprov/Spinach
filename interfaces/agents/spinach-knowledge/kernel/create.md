# kernel/create.m

Source: [kernel/create.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/create.m)
Wiki: [Spin Dynamics Wiki: create.m](https://spindynamics.org/wiki/index.php?title=create.m)

- Signature: `spin_system=create(sys,inter)`

## Purpose and inputs

The kernel entry point constructs the `spin_system` object used by the rest of Spinach. It validates and absorbs system and interaction specifications, then reports diagnostics.

- `sys` — spin-system and instrument specification structure. The source requires isotope labels and a scalar magnet field; consult the Spin System Specification section of the manual for the complete input specification.
- `inter` — interaction specification structure. It may be omitted; the executable code then sets it to `[]` before validation.
- `spin_system` — assembled system object, including system, component, interaction, chemistry, and relaxation data. For `N=numel(sys.isotopes)`, particle metadata is held per particle and pairwise couplings in `N-by-N` cells, with populated spin-coupling entries represented by `3x3` tensors.

## Assembly in execution order

1. Startup checks and configuration establish output and scratch destinations, run the available validation hooks, configure tolerances and enabled features, and prepare parallel/cache/GPU-related settings when requested.
2. Particle metadata is absorbed: isotope labels and particle count, spin versus bosonic-mode classification, multiplicities, gyromagnetic ratios, magnet field, and base frequencies `-gamma*magnet`.
3. Optional bosonic-mode specifications populate per-mode frequency, carrier, anharmonicity, damping, thermal occupation, and dephasing data, then pair-coupling and modulation data. If modes are present without mode interaction data, the source reports that zero couplings are assumed.
4. One-particle magnetic terms are assembled from susceptibility contributions and the supplied Zeeman eigenvalue/Euler-angle, matrix, or scalar forms; electron g-tensors and nuclear ppm shifts are converted using the base frequencies. Giant-spin coefficients are absorbed where present.
5. Chemistry, concentrations, exchange/flux settings, and order matrices are absorbed, with defaults where the source specifies them.
6. Coordinates and optional periodic-boundary vectors are absorbed and passed to `dipolar` for dipolar couplings. If coordinates are absent, the source assumes zero point-dipolar interactions and initialises an identity proximity matrix. User-specified pair couplings are then added.
7. Relaxation and radical-pair recombination specifications are absorbed and reported.
8. The final `inter.ignore` processing drops the listed coupling tensors, and leftover unparsed `sys` fields are reported as an error.

## Couplings, decoupling, and guards

`inter.ignore` is specifically a coupling drop list: for each listed pair the source clears both `{i,j}` and `{j,i}` entries after assembly. It does not remove spins or other interaction families. Pair couplings below the configured cutoff are cleared; couplings involving multiplicity-one ghost spins are also cleared. A populated coupling between different chemical subsystems is rejected.

Field validation is delegated substantially to `grumble` and the helper routines. The source checks the required isotope/magnet inputs and applies field-specific checks for modes, tensors, coordinates and periodic boundaries, chemistry, relaxation, recombination, and ignore-pair indices. In particular, damped bosonic modes require an explicit `inter.temperature`; bosonic particle combinations are restricted for mode couplings and modulation; and an unconsumed system option is an error. This function assembles the model specification; it does not solve an eigenfield problem or compute numerical derivatives.

## Units and conventions

The magnet field is in tesla; coordinates and periodic translations are in ångströms; temperatures are in kelvin; mode lifetimes, mode `T2` values, and correlation times are in seconds. Frequency-specified mode frequencies, carriers, anharmonicities, mode couplings and modulation terms, giant-spin coefficients, and user spin couplings are multiplied by `2*pi` on absorption. The coupling cutoff is applied to input coupling magnitudes in Hz. Mode lifetime damping is `1/lifetime`, linewidth damping is `2*pi*linewidth`, and Q-factor damping uses physical angular frequency divided by Q. Nuclear shifts are specified in ppm; electron Zeeman tensors use g-tensor units; absorbed Zeeman matrices are in angular-frequency units. The source reports relaxation and recombination rates in Hz. Consult the source-specific field validation for each interaction family's allowed forms.

For bosonic modes, `inter.modes.carriers` declares the laboratory rotating-frame frequency used by `inter.modes.frqs`; those declared frequencies are detunings and may be negative. Thermal occupation uses the physical frequency: carrier plus detuning when a positive carrier is specified, otherwise the absolute declared frequency. The mode pure-dephasing rate is `1/T2-kappa*(1+2*nbar)/2`, with amplitude-damping rate `kappa` and thermal occupation `nbar`; the source checks that the resulting pure-dephasing contribution is physically non-negative. Bosonic-mode quadratures use `(a+a')/sqrt(2)` in both `inter.modes.longitudinal` and the modulation channels.

## Call syntax

The source gives the call syntax `spin_system=create(sys,inter)`; when there are no interaction specifications, `create(sys)` is also accepted by the executable code.

Zero track elimination is off by default. Add `'zte'` to `sys.enable` to opt in; `'zte'` is no longer accepted in `sys.disable`. Paranoia overrides the opt-in and leaves ZTE off.
