# examples/nmr_spen/ufdosy_2spin.m

- Signature: `ufdosy_2spin()`

## Purpose

Simulates ultrafast DOSY for two coupled spins, including dipole-dipole (DD) and chemical-shift-anisotropy (CSA) relaxation, diffusion, and spatial flow. The source estimates minutes on an NVIDIA Tesla A100 and much longer on CPU. Authors: Ludmilla Guduff, Jean-Nicolas Dumez, and Ilya Kuprov.

## Model and setup

The system is two 1H spins at 14.1 T, with shifts 6.5 and 7.5, a 15 Hz scalar coupling, and CSA tensors with eigenvalues [-10 -10 20] and Euler angles [0 0 0] and [0 pi/2 0]. Redfield relaxation uses tau_c=1.0e-9 s, secular retention, and zero equilibrium. The basis is sphten-liouv with no approximation.

The one-dimensional sample model has length 0.015 m and 3000 points; flow velocity is 1e-4 and diffusion coefficient is 8e-10 m^2/s. The sequence calls imaging with spendosy; acquisition uses 128 points, 256 loops, deltat=1.5e-6 s, and Ga=0.52 T/m. Encoding uses 1000 pulse points, Te=0.0015 s, Tau=0.0016 s, BW=110000 Hz, and Ge=0.2535 T/m.

## Processing

The returned signal is Fourier transformed along both dimensions and plotted as magnitude against chemical shift and field of view.
