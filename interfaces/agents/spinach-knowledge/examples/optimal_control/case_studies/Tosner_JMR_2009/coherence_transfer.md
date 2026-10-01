# examples/optimal_control/case_studies/Tosner_JMR_2009/coherence_transfer.m

> This historical Spinach variant optimises a ten-level RF-power ensemble and uses L-BFGS; the paper’s first example is a single-system two-channel design. Its on-resonance 1H–13C state transfer and 1/J duration remain useful, but it does not reproduce the paper’s pulse.

- Signature: `coherence_transfer()`

## Purpose

Design a finite-duration heteronuclear coherence-transfer pulse from proton transverse magnetisation to carbon transverse magnetisation, Hx -> Cx. The source header identifies this as the first optimal-control example and cites the journal article by DOI; the function contains a pulse-design specification, not a reported optimisation result.

## Physical / mathematical content

The model is a two-spin 1H-13C system at 14.1 T. Both isotropic chemical shifts are set to 0 ppm, and the scalar coupling is 140 Hz. The NMR drift Hamiltonian is built from this system. The initial density operator is Lx on spin 1 (1H), and the target is Lx on spin 2 (13C); each is normalised. The fixed transfer time is T = 1/J, where J is the specified coupling.

## Numerical / algorithmic content

The complete `sphten-liouv` basis is selected with approximation `none`. Four Cartesian transverse controls are used: Lx and Ly on each isotope, mapped as two RF channels. The 1/J interval is divided into 150 equal slices, each of duration 1/(140*150) seconds. The available power levels are `2*pi*linspace(10,1000,10)` in angular-frequency units. The optimisation uses `fmaxnewton` with `@grape_xy`, the `lbfgs` method, at most 200 iterations, and an `NS` penalty of weight 0.01. Its initial guess is `rand(4,150)/10`; no random seed is set in the function. Requested diagnostics are XY controls, a spectrogram, and robustness plots.

## Syntax

Call `coherence_transfer()` with no arguments. Source: [examples/optimal_control/case_studies/Tosner_JMR_2009/coherence_transfer.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/case_studies/Tosner_JMR_2009/coherence_transfer.m).

## Parameters / inputs

The system, coupling, transfer interval, control channels, discretisation, power levels, penalty, iteration limit, and initial-guess dimensions are fixed in the function; there are no function inputs.

## Outputs

The `fmaxnewton` call is not assigned to an output variable. The function requests control, spectrogram, and robustness visualisations; it does not encode a numerical transfer fidelity or a returned pulse in its interface.

## Header notes

The source cites [DOI 10.1016/j.jmr.2008.11.020](https://doi.org/10.1016/j.jmr.2008.11.020). No article title or experimental result is asserted here.
