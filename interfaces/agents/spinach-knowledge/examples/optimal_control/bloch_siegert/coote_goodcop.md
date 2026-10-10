# examples/optimal_control/bloch_siegert/coote_goodcop.m

## Purpose

This example sets up a 13C optimal-control pulse design based on the GOODCOP logic cited by [Coote et al.](https://doi.org/10.1038/s41467-018-05400-4). The design applies a 150 us pulse that gives C-alpha spins their offset-dependent free evolution for 90% of the pulse duration while inverting longitudinal magnetisation across the CO band. The example links the duration, RF ceiling, contraction factor, and ensemble bands to Table 1 and the paper text, and the 53.2 ppm carrier to Supplementary Figure 5.

## Model and target

The model is one 13C spin at 18.8 T, corresponding to an 800 MHz 1H field, in the sphten-liouv basis with no basis approximation. The carrier is 53.2 ppm. The model uses normalised Ix, Iy, and Iz state vectors and the NMR drift Hamiltonian.

The C-alpha ensemble contains 100 uniformly spaced offsets from 35 to 75 ppm. Its initial transverse states alternate between Ix and Iy; each target is the state propagated by the offset Hamiltonian 2*pi*ca_hz*Lz for 0.90*150 us, or 135 us. The CO ensemble contains 30 uniformly spaced offsets from 165 to 185 ppm, with Iz as each initial state and -Iz as its target. Chemical-shift offsets are converted to Hz relative to the carrier at the specified field.

## Controls and optimisation

The pulse has two 13C controls, along x and y, sampled in 75 equal intervals of 2 us. The RF ceiling is 15 kHz; the control amplitude scale is 2*pi*15 kHz in angular-frequency units. The offset ensemble is associated with Lz, and the ensemble correlation is rho_ens. The control setup specifies the lbfgs method and a maximum of 200 iterations; fmaxnewton is called with the grape_xy objective.

Two waveforms are optimised from the same randomly drawn initial guess: one with Bloch-Siegert corrections enabled in optimcon, and one with those corrections disabled. The returned control arrays are multiplied by the RF amplitude scale. The Bloch-Siegert-enabled simulation augments the x/y controls with virtual controls using bloch_siegert; the uncorrected waveform is simulated without that augmentation.

## Evaluation

The evaluation grid contains 251 uniformly spaced chemical shifts from 0 to 200 ppm. At each point, the drift includes the corresponding Lz offset; Iz is propagated through each pulse, and the real final Iz overlap is recorded as Mz. The plotted observable is final Mz versus 13C chemical shift, with separate profiles for Bloch-Siegert corrections on and off.

## Source

[coote_goodcop.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/bloch_siegert/coote_goodcop.m)
