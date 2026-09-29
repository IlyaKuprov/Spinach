# examples/quantum_tech/jaynes_cummings_a.m

[Example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/jaynes_cummings_a.m)

- Signature: `jaynes_cummings_a()`

## Model and basis

This is a closed, coherent Jaynes–Cummings cavity-QED eigenvalue model: one electron spin-1/2 couples to a quantised cavity mode. The cavity is represented by the `C5` truncation (five population levels, as stated in the example header); no Tavis–Cummings ensemble is present. The script sets the magnetic field to 0.33 T, sets the cavity frequency from the electron Zeeman frequency so the two are resonant at zero detuning, and places the spin and cavity coupling in `inter.modes.exchange{1,2}=2.828e6` in Spinach's angular-frequency Hamiltonian convention. It switches to the cavity rotating frame with `assume(spin_system,'cavity')`.

The Hilbert-space basis is `zeeman-hilb` with no approximation. The electron detuning term is formed with its `Lz` operator. The calculation then projects onto the two one-excitation states: the excited-spin/empty-cavity state and the ground-spin/one-photon state. It diagonalises this two-state Hamiltonian at each detuning. There is no time-dependent drive or dissipative term in this source.

## Avoided crossing

The detuning sweep covers −15 to +15 MHz in 100 points (the source uses angular frequency internally and divides by `2*pi` for the plotted axis). The observable is the pair of one-excitation energy eigenvalues, plotted in MHz against detuning in MHz. Their avoided crossing is the model's vacuum-Rabi splitting from coherent spin-cavity coupling; it is a Hamiltonian spectrum, not a measured cavity trace or a device-fidelity result.
