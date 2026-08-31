# kernel/grids/ngridpts.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/ngridpts.m`
- Signature: `n=ngridpts(grad_amps,grad_durs,isotope,max_coh_order,sample_size)`
- Total lines: 83

## Purpose

Estimates the minimum number of spatial grid points necessary to have a valid treatment of gradient driven experiments with expli- cit digitization of spatial dimensions. Syntax: n=ngridpts(grad_amps,grad_durs,isotope,... max_coh_order,sample_size);

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 39-40: Check consistency; implemented by `grumble(grad_amps,grad_durs,isotope,max_coh_order,sample_size)`.
- Lines 42-43: Compute the worst-case total effective gradient; implemented by `worst_case_grad=sum(abs(grad_amps.*grad_durs))`.
- Lines 45-46: Compute the spatial frequency of the worst-case spiral; implemented by `worst_case_freq=abs(max_coh_order*spin(isotope)*worst_case_grad)`.
- Lines 48-49: Determine the minimum number of points required; implemented by `n=ceil(worst_case_freq*sample_size/pi)`.

### Key state/data transformations

- Lines 43: computes `worst_case_grad` using `worst_case_grad=sum(abs(grad_amps.*grad_durs))`.
- Lines 46: computes `worst_case_freq` using `worst_case_freq=abs(max_coh_order*spin(isotope)*worst_case_grad)`.
- Lines 49: computes `n` using `n=ceil(worst_case_freq*sample_size/pi)`.

### Local helper functions

- Line 54: `grumble()` — `function grumble(grad_amps,grad_durs,isotope,max_coh_order,sample_size)`.
  - Representative operation: `if (~isnumeric(grad_amps))||(~isreal(grad_amps))||(~isrow(grad_amps))`.
  - Representative operation: `error('grad_amps must be a row vector of real numbers.')`.

## Parameters / inputs

- grad_amps -a row vector of all gradient amplitudes
- in the sequence, T/m
- grad_durs -a row vector of all gradient durations
- in the sequence, s
- isotope -the highest magnetogyric ratio isotope in
- the spin system, e.g. '1H'
- max_coh_order -maximum order of coherence (either positive
- or negative) expected during the experiment
- being simulated
- sample_size -spatial extent of the sample, m

## Outputs

- n -the minimum recommended number of discretisation points
- Note: the function returns the minimum number of points, it may
- in practice be necessary to have several times the number,
- depending on your accuracy requirements.

## Implementation structure

- Estimates the minimum number of spatial grid points necessary to
- have a valid treatment of gradient driven experiments with expli-
- cit digitization of spatial dimensions. Syntax:
- n=ngridpts(grad_amps,grad_durs,isotope,...
- max_coh_order,sample_size);
- grad_amps -a row vector of all gradient amplitudes
- in the sequence, T/m
- grad_durs -a row vector of all gradient durations
- in the sequence, s
- isotope -the highest magnetogyric ratio isotope in
- the spin system, e.g. '1H'
- max_coh_order -maximum order of coherence (either positive

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `isrow()`, `any()`, `ischar()`, `isscalar()`.
