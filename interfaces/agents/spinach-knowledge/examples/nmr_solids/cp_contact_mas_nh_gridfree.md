# examples/nmr_solids/cp_contact_mas_nh_gridfree.m

- Signature: cp_contact_mas_nh_gridfree()

## Purpose

A simulated 1H-to-15N cross-polarisation (CP) contact curve for one 15N and one 1H, in the doubly rotating frame. This is the grid-free Fokker-Planck treatment of a spinning powder, initialised from thermal equilibrium. The source estimates minutes on a Tesla A100 GPU and much longer on a CPU.

## Spin model and experiment

The two-spin model has zero isotropic Zeeman shifts, coordinates [0, 0, 0] and [0, 0, 1.05] (the existing page identifies the separation as 1.05 angstrom), and temperature set to 298. The example specifies the shifts and coordinates rather than listing coupling tensors. It uses the full sphten-liouv basis with no approximation.

MAS is represented by a rotor-rate setting of 10000 and axis [sqrt(2/3), 0, sqrt(1/3)]; the source does not state the unit for that rate. The grid-free propagator uses max_rank 42. It requests iso_eq, detects the 15N Lx state, and supplies 100 time steps of 1e-5 seconds each. The two spin-lock nutation-frequency settings are 5e4 Hz on 1H and 4e4 Hz on 15N (50 and 40 kHz); these unequal values are the code inputs, not a separate assertion that a Hartmann-Hahn match was measured. The excitation operators are Hx on 1H and Ny on 15N; the spin-lock operators are Hy on 1H and Nx on 15N.

The example calls gridfree with the generic cp_contact_hard experiment function. That function applies the specified ideal pi/2 excitation and evolves under the supplied spin-lock terms; the example supplies the operators, powers, and time grid. No HMQC transfer/reconversion or DOR rotor sequence is configured here.

## Output and interpretation

The returned simulated signal is plotted as its real part against cumulative time in seconds, labelled as the 15N SX expectation value. It is not a measured spectrum. The source comment says a GPU is needed, but the sys.enable={'gpu'} line is commented out; the example itself does not actively set that option.

## Source

https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/cp_contact_mas_nh_gridfree.m
https://github.com/IlyaKuprov/Spinach/blob/main/experiments/cp_contact_hard.m
