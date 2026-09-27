# examples/optimal_control/mas_powder_control.m

- Signature: `mas_powder_control()`

## Purpose

Designs a phase-modulated pulse that transfers Lz to Ly on 87Rb in a quadrupolar rubidium system under magic angle spinning. Calculation time: hours.

## Physical / mathematical content

- The initial and target states are Lz and Ly, respectively. The ensemble includes powder orientations, three RF power levels, and five offsets from −1 to +1 kHz.
- The phase waveform is optimised using `fmaxnewton` with `grape_phase` and the limited-memory BFGS (`lbfgs`) method. The amplitude profile is fixed.

## Numerical / algorithmic content

- The optimised phase waveform is converted to Cartesian RF waveforms. A `parfor` loop propagates the pulse for each drift Liouvillian and computes the average target-state fidelity. The source does not specify GPU execution.

## Implementation structure

- Specify the 87Rb spin system, quadrupolar coupling, basis, and MAS parameters.
- Generate the ensemble drift Liouvillians and define the initial and target states.
- Define the Lx and Ly control operators, Lz offset operator, pulse timing, power levels, offsets, and initial phase guess.
- Optimise the pulse, evaluate its fidelity across the ensemble, and report the average fidelity.
