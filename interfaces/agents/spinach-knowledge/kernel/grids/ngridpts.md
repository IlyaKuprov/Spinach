# kernel/grids/ngridpts.m

- Signature: `n=ngridpts(grad_amps,grad_durs,isotope,max_coh_order,sample_size)`

## Purpose

Estimates the minimum number of spatial grid points for explicit spatial discretisation of gradient-driven experiments.

## Physical / mathematical content

The worst-case total effective gradient is `sum(abs(grad_amps.*grad_durs))`. The worst-case spatial frequency is `abs(max_coh_order*spin(isotope)*worst_case_grad)`.

## Numerical / algorithmic content

The function returns `n=ceil(worst_case_freq*sample_size/pi)`. It checks the inputs before calculating this value.

## Parameters / inputs

- `grad_amps`: row vector of gradient amplitudes in the sequence, T/m.
- `grad_durs`: row vector of gradient durations in the sequence, s; each duration must be positive, and the vector must have the same number of elements as `grad_amps`.
- `isotope`: character string naming the isotope with the highest magnetogyric ratio in the spin system, e.g. `'1H'`.
- `max_coh_order`: integer maximum coherence order, positive or negative, expected during the simulated experiment.
- `sample_size`: positive spatial extent of the sample, m.

## Outputs

- `n`: minimum recommended number of discretisation points. Several times this number may be needed, depending on accuracy requirements.

## Implementation structure

A local `grumble` function checks input types, shapes, positivity, and matching vector lengths. The main function then computes the worst-case gradient and spatial frequency and rounds the resulting point count up with `ceil`.

[Function documentation](https://spindynamics.org/wiki/index.php?title=ngridpts.m)