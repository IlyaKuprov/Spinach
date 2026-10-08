# examples/nmr_solids/mas_powder_nqi_floquet.m

- MATLAB implementation: [examples/nmr_solids/mas_powder_nqi_floquet.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_nqi_floquet.m)

[MATLAB source](../../../../../examples/nmr_solids/mas_powder_nqi_floquet.m)

## Purpose and spin system

This example computes a powder MAS spectrum for one quadrupolar `2H` nucleus. It sets the field to `9.4 T`, the quadrupolar tensor eigenvalues to `[-1e3 -2e3 3e3] Hz`, and its Euler angles to `[0 0 0]`. The values are model inputs, not experimental measurements; the tensor eigenvalues sum to zero. Spinach documents quadrupolar interaction tensors in Hz in its [g2spinach knowledge page](../../interfaces/g2spinach.md).

## MAS algorithm and acquisition

The source uses the Floquet route directly: `floquet(spin_system,@acquire,parameters,'nmr')`. Its header says perturbative corrections to the rotating-frame transformation are not applied and estimates a seconds-scale calculation time; the timing is not a measured runtime. The rotor axis is `[1 1 1]` and rate `1000 Hz`; the powder grid is `leb_2ang_rank_17` with `max_rank=17`. Acquisition uses a `2e4 Hz` sweep, 512 points, zero-fill 4096, and offset 0. It selects `2H`, leaves `decouple={}`, sets ppm axis units, and inverts the axis. Both initial state and receiver are `L+` on `2H`.

The returned FID is exponentially apodised with parameter `6`, Fourier transformed using `fftshift(fft(fid,parameters.zerofill))`, and the real spectrum is plotted with `plot_1d`. This is a computed spectrum; the source does not provide experimental measured output or a numerical comparison.
