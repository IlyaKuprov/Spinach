# examples/relaxation_theory/trosy_nh.m

- Source: [examples/relaxation_theory/trosy_nh.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/trosy_nh.m)
- Signature: `trosy_nh()`

## Purpose

Calculate transverse relaxation matrix elements across magnetic field for a two-spin amide `1H–15N` model. The source states a 25 ns rotational correlation time, nitrogen CSA parameters from [the cited study](https://doi.org/10.1021/ja0016194), and an N–H bond length taken from DFT. It estimates a runtime of minutes.

## Model and relaxation pathway

The spin system contains `1H` and `15N`. The source supplies Zeeman eigenvalues [6, 0, −6] for the proton and [−108, 62, 46] for nitrogen, with Euler-angle inputs [0, 0, 0] and [0, 0, −19] degrees (converted in the code to radians). Coordinates are [1.04, 0, 0] and [0, 0, 0]; the source does not state a unit for these coordinates. The tensors and internuclear geometry define the anisotropic-shielding and dipolar contributions used by the relaxation calculation; the source does not print a component-by-component decomposition of the resulting superoperator.

The script explicitly calls the Redfield relaxation model (`inter.relaxation={'redfield'}`), keeps relaxation in the lab frame, sets equilibrium to zero, and uses `tau_c={25e-9}` (25 ns). It uses the `sphten-liouv` formalism with no basis approximation. For the TROSY-style comparison it evaluates normalised single-spin transverse coherences and paired operators with opposite-sign two-spin terms: proton raising coherence with proton-plus-nitrogen-longitudinal coherence, and nitrogen raising coherence with nitrogen-plus-proton-longitudinal coherence. These left/right branches expose the relaxation interference relevant to the TROSY comparison; the script does not separately label or export an individual cross-correlation term. It does not call a stochastic-Liouville solver.

## Calculation and output

The script samples 30 proton Larmor frequencies from 200 to 1500 MHz, converts each to a field with `2*pi*lin_freq*1e6/spin('1H')`, rebuilds the spin system and basis, and evaluates the relaxation superoperator at that field. It projects the superoperator onto the single-spin and paired operators; the plotted vertical quantity is a relaxation matrix element in Hz. Two figures show the proton and nitrogen operator matrix elements, respectively; each compares the single-spin transverse coherence with its two opposite-sign branches against proton Larmor frequency. This is a calculation of operator relaxation elements, not a pulse sequence, acquired signal, or measured spectrum.
