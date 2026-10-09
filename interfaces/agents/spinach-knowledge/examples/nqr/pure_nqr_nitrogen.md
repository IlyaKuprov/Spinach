# examples/nqr/pure_nqr_nitrogen.m

- Signature: `pure_nqr_nitrogen()`
- Source: [examples/nqr/pure_nqr_nitrogen.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nqr/pure_nqr_nitrogen.m)

## Objective and model

This example specifies a zero-field powder NQR calculation for one spin-1 `14N` nucleus. It sets `sys.magnet=0` and builds the quadrupolar interaction with `eeqq2nqi(1.18e6,0.53,1,[0 0 0])`: the helper defines the first argument as `e^2 q Q / h` in Hz, the second as the dimensionless asymmetry, the third as the nuclear spin, and the last as Euler angles in radians. Thus the input is 1.18 MHz, asymmetry 0.53, spin 1, and zero orientation angles. The helper's convention is documented in [kernel/conventions/transforms/eeqq2nqi.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/eeqq2nqi.m).

The calculation uses the exact `sphten-liouv` basis with `bas.approximation='none'`. Relaxation is configured with `inter.relaxation={'damp'}`, `inter.damp_rate=1e5`, lab-frame retention, zero equilibrium, and temperature 298 K.

## Acquisition and processing

The acquisition requests anisotropic equilibrium, uses a 5e6 Hz sweep, 512 points, and the `rep_2ang_200pts_sph` powder grid. The receiver operator is the `14N` `L+` state; the pulse is an `Lx` operator with flip angle `pi/2`. Spinach's `powder` acquisition calls `hp_acquire` in the lab frame. The returned FID is exponentially apodised with parameter 6, then transformed as `imag(fftshift(fft(fid)))`; plotting uses MHz axis units.

Estimated calculation time: seconds.
