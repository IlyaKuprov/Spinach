# examples/nmr_solids/hmqc_mas_dq.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/hmqc_mas_dq.m

## What it models

The source describes rotor-synchronised powder MAS CN2D detection for a 14N–1H pair, using a one-dimensional Fokker–Planck treatment on a spherical orientation grid. Its stated target includes the second-order quadrupolar shift and lineshape. The Hamiltonian setup assigns the field as 19.96, an 14N quadrupolar interaction through eeqq2nqi(1.18e6, 0.53, 1, [0 0 0]), and Zeeman scalars 32.4 and 5 for the two spins. The source sets the coordinates to [0 0 0] and [0 0 1]. The earlier page rendered the field and quadrupole value as 19.96 T and 1.18 MHz; the source's numeric assignments are retained here, while its code does not spell out units alongside them. No explicit 14N–1H coupling matrix is assigned in this file.

The basis is sphten-liouv with approximation none. The experiment setup assigns rate 125000 about [1 1 1], max_rank 16, and grid rep_2ang_200pts_sph. It sets sweep [rate/4, 20000], 256 by 128 acquired points and zero-fill [1024, 512], uses 14N and 1H with the 14N rotating frame at rank 3, and labels the axes ppm. The initial state is a transverse 1H combination and the receiver is the 1H L+ operator. RF settings are rf_pwr=40e3 and rf_dur=2e-3; these assignments have no units stated in the source. There is no explicit Hartmann–Hahn match or contact-time assignment in this wrapper.

## Sequence wrapper and output

The wrapper calls singlerot with @cn2d_dq and the qnmr mode. That selects the DQ sequence callback but does not define its pulse/transfer internals in this example; those details should not be inferred from the wrapper. The cosine and sine FIDs are each apodised with sqcos in both dimensions, Fourier transformed in the two dimensions, combined as a States signal, and the real spectrum is passed to plot_2d.

The source comment estimates hours on CPU or minutes on a Tesla V100 GPU. This is a code comment, not a measured result reported by this page; the example does not supply a measured spectrum, benchmark record, or validation claim.
