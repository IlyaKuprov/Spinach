# examples/nmr_liquids/ct_cosy_rotenone.m

- MATLAB implementation: [examples/nmr_liquids/ct_cosy_rotenone.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/ct_cosy_rotenone.m)

- Signature: `ct_cosy_rotenone()`

## Purpose

A liquid-state constant-time COSY simulation for the 22-proton rotenone assignment cited in [DOI 10.1002/jhet.5570250160](http://dx.doi.org/10.1002/jhet.5570250160). The source estimates a calculation time of minutes.

## Spin system and sequence

All 22 sites are 1H. Their ordered chemical shifts are [6.72, 6.40, 4.13, 4.56, 4.89, 6.46, 7.79, 3.79, 2.91, 3.27, 5.19, 4.89, 5.03, 1.72, 1.72, 1.72, 3.72, 3.72, 3.72, 3.76, 3.76, 3.76] ppm at 5.9 T. The listed scalar couplings are sparse rather than all-to-all: examples include J(3,4)=12.1 Hz, J(9,10)=15.8 Hz, J(10,11)=9.8 Hz, and J(7,9)=J(7,10)=0.7 Hz; the source also specifies the other listed pair couplings. The assignment is a model input, not a fitted spectrum in this script.

The basis uses the `sphten-liouv` formalism, the `IK-2` approximation with scalar-coupling connectivity and proximity level 1, greedy system reduction, and a proximity cutoff of 4.0. Three S3 symmetry groups are assigned to sites [14,15,16], [17,18,19], and [20,21,22]. The sequence call is `liquid(...,@ct_cosy,...,'nmr')`; it selects 1H, sets the pulse-angle parameter to pi/2, offset to 1200, sweeps to [2000 2000], 256 points and 512 zero-fill points per dimension, and requests ppm axes.

No separate mixing-time, phase-cycle table, or receiver-phase setting appears in this wrapper. Those details are not independently established by this example; they belong to the called `ct_cosy` implementation and its conventions.

## Processing and output

A cosine window is applied along both dimensions, followed by a shifted 2D FFT. The plotted quantity is the spectrum magnitude in positive mode. This is a calculated spectrum for the cited assigned spin system; the script does not load a measured rotenone spectrum or provide a simulated-versus-experimental comparison. Its reduced-basis and symmetry settings are part of this example's computational model.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
