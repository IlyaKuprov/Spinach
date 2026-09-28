# examples/quantum_tech/spin_phonon_swap.m

- Signature: `spin_phonon_swap()`

## Purpose

Models resonant excitation exchange between an electron spin and a quantised phonon mode in the spin–phonon Jaynes–Cummings limit used in mechanical spin-qubit proposals such as Rabl et al., Nature Physics 6, 602 (2010). Calculation time: seconds.

## Physical / mathematical content

- A single spin excitation is exchanged with the phonon mode. The example tracks both excitation populations, checks that transfer is visible, and verifies population conservation in the active doublet.

## Numerical / algorithmic content

- The `zeeman-hilb` model uses no basis approximation. The `cavity` device context preserves spin–mode exchange, and the trajectory contains 801 points over 0–800 ns.

## Implementation structure

- The system uses isotopes `{'E','V3'}`, a resonant phonon mode at zero rotating-frame frequency, and exchange coupling `4e6`. It starts in `{'ZL2','BL1'}` (spin excitation, phonon vacuum); projectors `{'ZL2','E'}` and `{'ZL1','BL2'}` measure spin and phonon populations.
