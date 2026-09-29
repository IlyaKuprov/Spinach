# examples/nmr_solids/mas_powder_dip_fplanck.m

Source: [examples/nmr_solids/mas_powder_dip_fplanck.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_dip_fplanck.m)

## Model

This is a simulated MAS powder pulse-acquire example for two dipole-coupled protons. The source comment cites DOI [10.1016/j.jmr.2016.07.005](https://doi.org/10.1016/j.jmr.2016.07.005) for its Fokker-Planck formalism. The system is set to 14.1 T; its `1H` isotropic Zeeman scalar values are 5.0 and -2.0, with positions [0, 0, 0] and [0, 3.9, 0.1]. The coordinates supply the dipolar geometry when `create` builds the system; there is no separate dipolar-coupling constant in the input. This file does not state units for the coordinates or Zeeman scalars. It does not include a quadrupolar spin or NQI, nor an explicit RF pulse or decoupling field.

The basis is spherical-tensor Liouville space with no approximation and projection +1. The MAS axis is [1, 1, 1] and the rotor rate is 1000 Hz; maximum rank is 17 and the configured grid is `leb_2ang_rank_17`. Unlike the companion Floquet example, this script calls `singlerot(spin_system,@acquire,parameters,'nmr')`, the source's Fokker-Planck pulse-acquire path.

## Signal and display

The initial state and coil are the proton `L+` state. The acquisition has 512 points, a sweep setting of 2e4, and zero filling to 4096. The resulting FID receives exponential apodisation parameter 6, is Fourier transformed with `fftshift`, and its real spectrum is plotted by `plot_1d`. This is the calculated spectrum from the model settings; the source does not present it as an experimental measurement.