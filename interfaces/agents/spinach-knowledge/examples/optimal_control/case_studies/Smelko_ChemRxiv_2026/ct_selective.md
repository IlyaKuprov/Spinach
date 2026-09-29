# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/ct_selective.m

## Objective and context
Design a central-transition-selective pulse for the z-filtered 27Al MQMAS experiment. Related ChemRxiv record: [DOI 10.26434/chemrxiv.15008427](https://doi.org/10.26434/chemrxiv.15008427).

## Spin model and rotor ensemble
The model is one spin-5/2 27Al nucleus with quadrupolar coupling CQ = 3.0 MHz, asymmetry eta = 1.0, and 10 ppm axial shielding anisotropy, represented by principal values [-5, -5, 10] ppm. It uses a 400 MHz proton-frequency reference, 12.5 kHz magic-angle spinning, and a second-order quadrupolar interaction in the rotating frame. Rotor-phase-resolved drifts span 200 crystallite orientations and 80 initial rotor phases, with 160 rotor ticks.

## State transfer and pulse design
The normalised initial operator is the population difference across the central transition. The target is a single-quantum coherence between the central-transition levels; the GRAPE objective is transfer fidelity to this target across the powder and rotor-phase drift ensemble. The pulse uses Cartesian Lx and Ly controls in 100 slices of 0.5 us (50 us total), with a 10 kHz amplitude ceiling. An SNSA amplitude-spillout penalty of weight 100 constrains an L-BFGS GRAPE search for up to 500 iterations. The initial guess has random amplitudes up to 10% of the ceiling with one slice at the ceiling. After optimisation, amplitudes are clipped to the ceiling and fidelity is reevaluated.

Reference-waveform fidelity after 500 iterations: 0.69; the maximum for this normalised initial/target state pair is 1/sqrt(2). The waveform is stored in rad/s with its slice durations in ct_pulse.mat.

## Source
[ct_selective.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/ct_selective.m)
