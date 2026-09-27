# examples/relaxation_theory/scalar_relaxation_1.m

- Signature: `scalar_relaxation_1()`

## Purpose

Constructs and displays a Redfield superoperator for scalar relaxation of the first kind in a two-proton system with a fluctuating J-coupling. The example is described as modelling the effect of slow nitrogen inversion in aziridines, which modulates scalar couplings on a millisecond timescale; see the cited [article](http://dx.doi.org/10.1002/ange.201410271). Calculation time: seconds.

## Physical / mathematical content

The two `1H` spins are specified at `11.75 T`, with Zeeman scalars `0.0` and `2.0`. The Redfield setup selects `SRFK` relaxation, keeps the `kite` terms, uses zero equilibrium, sets correlation times to `1.0` and `1e-3`, and sets the sole off-diagonal SRFK modulation depth to `15.0`.

## Numerical / algorithmic content

The model uses the complete `sphten-liouv` basis with no approximation. It creates the spin system and basis, obtains the relaxation superoperator, and displays the nonzero pattern with `spy`.
