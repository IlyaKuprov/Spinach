# examples/quantum_tech/jaynes_cummings_c.m

- Signature: `jaynes_cummings_c()`
- Source: [`examples/quantum_tech/jaynes_cummings_c.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/jaynes_cummings_c.m)

## Purpose

A time-domain two-electron, one-cavity calculation: each electron has a Jaynes–Cummings-type exchange with the shared mode, and the electrons also have a direct scalar exchange coupling. The source propagates initial transverse spin coherence with the cavity empty and plots the summed spin signal and a cavity-field quadrature.

## Physical model and spin selection

The declared system is `sys.isotopes={'E','E','C5'}` at `sys.magnet=0.33` T. Here `E` is Spinach's generic electron-spin isotope, while `C5` is a five-level cavity oscillator. The electron–electron scalar coupling is `5e6` Hz. The cavity frequency is set from the electron resonance expression `-sys.magnet*spin('E')/(2*pi)`; its two unequal spin–cavity exchange values are `2.828e6` and `2.728e6` Hz. Thus this is a shared-mode two-emitter Jaynes–Cummings extension (Tavis–Cummings-like), with unequal couplings and an additional direct spin exchange, rather than the uncoupled ideal Tavis–Cummings limit.

The source uses the `sphten-liouv` formalism with no basis approximation and calls `device(...,'cavity')`. It selects `parameters.spins={'E'}`, with offset `5e6` Hz, sweep `1e8` Hz, and `251` points. This selects modeled electron-spin transitions; no named defect, defect-specific nuclear isotope, hyperfine interaction, anisotropic g tensor, powder/orientation average, or measured EPR spectrum is specified.

## Initial state and plotted observables

The initial state is the sum of the two electrons' `Lx` coherences, each paired with cavity state `BL1` (the empty-mode level). Detection adds the two electron `Lx` operators and separately projects the cavity onto `(C-A)/2i`. The source plots the real parts against a `251`-sample time axis from 0 to 2.5 μs.

These curves are simulated spin-coherence and cavity-quadrature trajectories for the declared coupled model. The source sets no independent drive-amplitude parameter, cavity linewidth, or relaxation term; its cavity-device context comes through `device(...,'cavity')` and the listed sequence parameters. The plot is not a measured spectrum or device-fidelity assessment, and no quantitative vacuum-Rabi splitting is asserted here.
