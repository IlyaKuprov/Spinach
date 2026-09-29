# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mq_excitation.m

## Objective and context
Design the multiple-quantum excitation step for the z-filtered 27Al MQMAS experiment: convert longitudinal magnetisation Iz into a 3Q or 5Q coherence. The function defaults to 5Q and accepts mq_order = 3 or 5. Related ChemRxiv record: [DOI 10.26434/chemrxiv.15008427](https://doi.org/10.26434/chemrxiv.15008427).

## Spin model and rotor ensemble
The model is one spin-5/2 27Al nucleus with quadrupolar coupling CQ = 3.0 MHz, asymmetry eta = 1.0, and 10 ppm axial shielding anisotropy, represented by principal values [-5, -5, 10] ppm. It uses a 400 MHz proton-frequency reference, 12.5 kHz magic-angle spinning, and a second-order quadrupolar interaction in the rotating frame. Rotor-phase-resolved drifts span 200 crystallite orientations and 80 initial rotor phases, with 160 rotor ticks.

## State transfer and pulse design
The initial state is normalised Iz. The target is the normalised Hermitian combination of the +MQ and -MQ coherences between m = +mq_order/2 and m = -mq_order/2: the m = +/-3/2 pair for 3Q or the m = +/-5/2 pair for 5Q. The GRAPE objective is transfer fidelity to this target across the powder and rotor-phase drift ensemble. Cartesian Lx and Ly controls are optimised in 480 slices of 0.5 us (240 us total, three rotor periods), with a 100 kHz amplitude ceiling. An SNSA amplitude-spillout penalty of weight 100 constrains an L-BFGS GRAPE search for up to 500 iterations. The initial guess has random amplitudes up to 10% of the ceiling with one slice at the ceiling; the optimised amplitudes are clipped and fidelity is reevaluated.

Reference-waveform fidelities after 500 iterations: 0.72 (3Q) and 0.65 (5Q). The saved waveform is in rad/s with slice durations, in mq_exc_<order>q.mat.

## Source
[mq_excitation.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mq_excitation.m)
