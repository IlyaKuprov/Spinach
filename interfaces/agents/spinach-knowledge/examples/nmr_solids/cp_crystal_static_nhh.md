# examples/nmr_solids/cp_crystal_static_nhh.m

- Signature: cp_crystal_static_nhh()

## Purpose

A simulated static single-crystal 1H-to-15N cross-polarisation (CP) contact curve for one 15N and an eight-proton bath. The source says it retains the full Liouville space because this is not a powder calculation and every spin interacts with every other spin. It estimates minutes on a Tesla A100 and longer on a CPU.

## Spin model and experiment

The nine-spin model specifies isotropic Zeeman shifts and explicit coordinates; the source describes the eight protons as scattered on a 2 angstrom-radius sphere around 15N. Temperature is set to 298. The basis is sphten-liouv with no approximation. The script specifies shifts and coordinates rather than listing coupling tensors.

This is one static crystal orientation, [pi/3, pi/4, pi/5], with no powder average or rotor. It requests aniso_eq for 15N, detects the 15N Lx state, and supplies 5e4 Hz (50 kHz) spin-lock nutation frequency on each channel. There are 100 time steps of 1e-5 seconds each. Excitation uses Hx on 1H and Ly on 15N; spin-lock uses Hy on 1H and Lx on 15N.

The call to crystal uses the generic cp_contact_hard experiment function: the helper applies the supplied ideal pi/2 excitation and evolves under the supplied spin-lock terms to return a contact curve. The example—not the helper—sets this spin model, full basis, orientation, RF values, and time grid. No HMQC transfer/reconversion or DOR sequence is configured.

## Output and interpretation

The returned simulated signal is plotted as its real part against cumulative time in seconds, labelled as the 15N SX expectation value. This is simulated output, not a measured spectrum. Although a source comment says a GPU is needed, the active setting is sys.enable={'greedy'}; 'gpu' appears only in a comment.

## Source

https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/cp_crystal_static_nhh.m
https://github.com/IlyaKuprov/Spinach/blob/main/experiments/cp_contact_hard.m
