# kernel/grids/ngridpts.m

- Signature: `n=ngridpts(grad_amps,grad_durs,isotope,max_coh_order,sample_size)`

## Purpose

Estimates the minimum number of spatial grid points necessary to have a valid treatment of gradient driven experiments with expli- cit digitization of spatial dimensions. Syntax: n=ngridpts(grad_amps,grad_durs,isotope,... max_coh_order,sample_size);

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

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
