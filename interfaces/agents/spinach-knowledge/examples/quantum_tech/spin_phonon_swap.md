# examples/quantum_tech/spin_phonon_swap.m

- Signature: `spin_phonon_swap()`
- Source: [MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/quantum_tech/spin_phonon_swap.m)

## Model

This is a resonant spin–phonon Jaynes–Cummings-limit excitation-swap model of the type discussed in mechanical spin-qubit proposals, including Rabl et al., *Nature Physics* **6**, 602 (2010). The particle list `{'E','V3'}` is a multiplicity-2 electron (spin one-half) and a three-level phonon mode; `V3` is mode truncation, not a vanadium or defect isotope. With `sys.magnet=0`, the phonon rotating-frame frequency is zero and the exchange coefficient is `4e6` (4 MHz). No particular defect nuclear isotope, higher-spin defect Hamiltonian, EPR field/frequency selection, external drive, or dissipation is specified.

## Initial state and plotted observable

Using `zeeman-hilb` with no basis approximation, the initial state `{'ZL2','BL1'}` places one excitation on the electron and none in the phonon mode. The `cavity` device context retains the spin–mode exchange. The source propagates 801 points with `sweep=1e9` and plots the projected real spin and phonon populations, measured with `{'ZL2','E'}` and `{'ZL1','BL2'}`, against an explicit 0–800 ns time axis. It checks for visible population transfer and conservation of the combined population in the active doublet.

## Interpretation

The curves represent ideal resonant exchange between one modelled electronic excitation and one phonon quantum. They are neither an EPR spectrum nor a measured spin-qubit swap or fidelity result. The source's “Calculation time: seconds” note is not measured here and does not establish runtime or convergence.
