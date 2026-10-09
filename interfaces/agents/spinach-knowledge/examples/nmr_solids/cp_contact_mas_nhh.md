# examples/nmr_solids/cp_contact_mas_nhh.m

- Signature: cp_contact_mas_nhh()

## Purpose

A simulated 1H-to-15N cross-polarisation (CP) contact curve for one 15N coupled to an eight-proton bath, in the doubly rotating frame. The example treats a spinning powder with a restricted Liouville-space basis retaining correlations through three spins. Its source estimates minutes on a Tesla A100 and much longer on a CPU.

## Spin model and experiment

The nine-spin model specifies the isotropic Zeeman shifts and explicit coordinates for the eight protons around 15N; the source describes the protons as scattered on a 2 angstrom-radius sphere. Temperature is set to 298. No explicit coupling tensors are listed in this example. It uses sphten-liouv with approximation IK-0 and inter_level 3. Interaction and proximity cutoffs are set to 5.0 and 4.0; the source gives no units for these cutoff values. It disables trajlevel and enables greedy.

The spinning-powder treatment uses the rep_2ang_100pts_sph grid and singlerot, with parameters.rate set to 10000; the source does not state the rate unit. max_rank is 3. The run requests iso_eq, observes the 15N Lx state, and uses 100 time steps of 1e-5 seconds each. Spin-lock nutation-frequency inputs are 5e4 Hz on 1H and 4e4 Hz on 15N (50 and 40 kHz); these are the configured values, not a claim that the two channels are Hartmann-Hahn matched. The code defines the spin-lock and excitation operators separately for the proton and nitrogen channels.

The example invokes singlerot with the generic cp_contact_hard experiment function. That function applies the supplied ideal pi/2 excitation and evolves under the supplied spin-lock terms; it does not supply this example's spin model, rotor/grid settings, RF values, or contact time. No HMQC transfer/reconversion or DOR sequence is configured.

## Output and interpretation

The returned simulated signal is plotted as its real part against cumulative time in seconds, labelled as the 15N SX expectation value. The example reports neither measured spectra nor a separate validation result.

## Source

https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/cp_contact_mas_nhh.m
https://github.com/IlyaKuprov/Spinach/blob/main/experiments/cp_contact_hard.m
