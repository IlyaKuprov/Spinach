# examples/quantum_tech/spin_cavity_vacuum_rabi.m

- Signature: `spin_cavity_vacuum_rabi()`
- Source: [MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/spin_cavity_vacuum_rabi.m)

## Model

This is a one-electron, one-cavity-mode Jaynes–Cummings calculation of vacuum-Rabi exchange, described in the source as the one-spin limit of the spin-ensemble cavity experiments of Schuster et al. and Kubo et al., *Physical Review Letters* **105**, 140501 and 140502 (2010). It is a schematic coherent exchange model, not a device-fidelity calculation or a measured defect spectrum.

The particle list `{'E','C3'}` means a multiplicity-2 electron (spin one-half) and a cavity mode truncated to three levels; `C3` is a mode specification, not a carbon isotope. The source sets `sys.magnet=0`, places the cavity at zero rotating-frame frequency, and sets its exchange coefficient to `8e6` (8 MHz in the frequency convention used by the source). The exchange-only model has no explicit drive, relaxation, or dissipation term. It is not a material-specific defect Hamiltonian: no nuclear defect isotope, defect zero-field splitting, or EPR field/frequency selection is specified.

## Initial state and plotted observable

In the `zeeman-hilb` formalism with no basis approximation, the initial state `{'ZL2','BL1'}` places one excitation on the electron and the cavity in its vacuum state. Spin and cavity excitation populations are projected from each trajectory point using `{'ZL2','E'}` and `{'ZL1','BL2'}`. The `cavity` device trajectory has 501 points and the figure plots the real populations over 0–500 ns. The source checks that transfer is visible and that the two populations sum to one in the active doublet; this is an internal numerical check, not experimental validation.

## Interpretation

The figure shows the ideal time-domain exchange of a single excitation between the two modeled subsystems. It does not show a driven EPR spectrum, an isotope-resolved defect transition, or loss-limited cavity dynamics. The source labels the calculation time as seconds, but that comment is not a measured runtime or a convergence claim.
