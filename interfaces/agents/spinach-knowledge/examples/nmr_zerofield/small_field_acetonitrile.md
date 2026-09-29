# examples/nmr_zerofield/small_field_acetonitrile.m

## Experiment represented

This is a simulated small-field 1H NMR spectrum for acetonitrile with 13C on the methyl carbon. The source says it is set to reproduce Figure 3 of [the cited Physical Review Letters paper](https://doi.org/10.1103/PhysRevLett.107.107601). The calculation generates a signal from the stated spin model; it does not load measured data.

## Spin model and field

The model contains three 1H spins and one 13C spin. The source labels the field as 2.64 mG and sets sys.magnet to 2.64e-3 × 1e-4 T, i.e. 2.64e-7 T. Each proton–carbon scalar coupling is 136.200 Hz. Temperature is 298 K. The basis is zeeman-hilb with approximation none. This nonzero small-field setting is distinct from the companion zero-field calculation; the code does not define a gradient or chirp.

## Acquisition and processing

The simulation uses a 700 Hz sweep, 4096 acquired points, and 16384 zero-fill points; offset is zero, the detected channel is 1H, axis units are Hz, and axis inversion is disabled. The nominal flip angle is π/2 and detection is uniaxial. Spinach builds the system and basis, then liquid propagates it with zerofield in the lab frame. The returned FID is mean-subtracted and exponentially apodised with parameter 6; the script computes the shifted FFT at the zero-fill length and plots its real part. The source header estimates a calculation time of seconds.

## Source

[MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_zerofield/small_field_acetonitrile.m)