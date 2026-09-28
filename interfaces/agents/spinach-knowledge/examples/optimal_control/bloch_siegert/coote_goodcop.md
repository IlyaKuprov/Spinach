# examples/optimal_control/bloch_siegert/coote_goodcop.m

- Signature: `coote_goodcop()`

## Purpose

Reproduces the GOODCOP pulse-design logic from [Coote et al.](https://doi.org/10.1038/s41467-018-05400-4), with Bloch-Siegert corrections enabled during design and validation. The pulse enforces contracted-time C-alpha evolution while inverting CO. Calculation time: minutes.

## Physical / mathematical content

- The single-(^{13}mathrm{C})-spin model is set at 18.8 T with a 53.2 ppm carrier. C-alpha offsets span 35–75 ppm (100 points) and CO offsets span 165–185 ppm (30 points). Alternating C-alpha initial states (I_x) and (I_y) are targeted to their drift-evolved states over the contracted duration; CO (I_z) states are targeted to (-I_z).
- The source identifies the duration, RF ceiling, contraction factor, two ensemble bands, and carrier as values from Table 1 and the paper text or Supplementary Figure 5. In the implementation these are `pulse_dur=150e-6` s, `max_rf_hz=15e3`, `alpha_scale=0.90`, and `carrier_ppm=53.2`.

## Numerical / algorithmic content

- The 150 μs pulse is represented by 75 slices. L-BFGS is limited to 200 iterations and uses the `rho_ens` ensemble correlation. The code designs a BSS-enabled pulse and a comparison pulse with BSS disabled, then evaluates both with BSS physics over a 251-point 0–200 ppm grid and plots final (M_Z).

## Implementation structure

- Set the single-spin model and sphten-liouv basis; construct control operators, normalized states, and drift Hamiltonian; define the paper parameters and offset ensembles; assemble correlated targets; optimize with BSS on and off; evaluate both waveforms on the dense chemical-shift grid and plot the GOODCOP inversion profile.
