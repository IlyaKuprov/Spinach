# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_efficiency.m

- Signature: `mqmas_efficiency()`
- Source: [examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_efficiency.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_efficiency.m)

## Purpose

This script compares simulated sequence efficiency for hard-pulse and optimal-control implementations of a z-filtered 27Al multiple-quantum magic-angle-spinning (MQMAS) experiment. It evaluates pulse waveforms made by the sibling excitation, conversion, and central-transition-selective optimisation examples; it does not perform those waveform optimisations itself. The source identifies the comparison with the ChemRxiv DOI [10.26434/chemrxiv.15008427](https://doi.org/10.26434/chemrxiv.15008427).

## Spin system and sequence

The model is one 27Al nucleus with spin 5/2, quadrupolar coupling 3.2 MHz, asymmetry 0.16, and a second-order quadrupolar interaction in the rotating frame. It includes shielding tensor eigenvalues -5, -5, and 10 ppm. The magnet is 400 MHz and the sample spins at 12.5 kHz about the axis [sqrt(2/3), 0, sqrt(1/3)]. The calculation uses a lab-frame Zeeman-Hilbert-space model.

The simulated order is selected as 3Q or 5Q. The sequence applies excitation, filters to the positive and negative selected multiple-quantum orders, applies conversion, filters to zero quantum order, applies a central-transition-selective pulse, and reads the central-transition single-quantum component. The initial state is normalised 27Al longitudinal magnetisation. Rotor-resolved drift Hamiltonians are sampled at 0.5 microsecond ticks; the powder average uses 400 crystallite orientations and 32 initial rotor phases per orientation.

## Pulse comparison and observable

The hard-pulse durations for excitation, conversion, and central-transition selection are [4.2, 1.4, 9] microseconds for 3Q and [4.4, 2.4, 9] microseconds for 5Q. Their angular-frequency amplitudes are 2*pi*[100, 100, 9.3] kHz. For the optimal-control comparison, the script loads the matching order-specific excitation and conversion waveforms and the central-transition waveform from the MAT files written by the sibling examples.

For each powder/rotor instance, the script propagates the density matrix through the three pulses and intervening coherence filters. It records matrix element (3,4), the central-transition single-quantum component, then defines efficiency as the modulus of the mean complex signal over the ensemble. Thus the reported quantity is not the mean of individual signal magnitudes.

## Source-stated comparison values

The source comments report 3Q efficiencies of 0.071 for hard pulses and 0.419 for optimal-control pulses, an enhancement factor of 5.9, compared with 5.7 simulated in the cited paper. For 5Q, the corresponding efficiencies are 0.0109 and 0.258, an enhancement factor of 23.8, compared with 25 simulated in the paper. The source comments estimate minutes of calculation time on 128 cores.
