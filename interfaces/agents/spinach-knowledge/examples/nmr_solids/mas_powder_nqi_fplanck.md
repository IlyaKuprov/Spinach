# examples/nmr_solids/mas_powder_nqi_fplanck.m

- MATLAB implementation: [examples/nmr_solids/mas_powder_nqi_fplanck.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_nqi_fplanck.m)

[MATLAB source](../../../../../examples/nmr_solids/mas_powder_nqi_fplanck.m)

## Purpose and spin system

This example computes a powder MAS spectrum for one quadrupolar `2H` nucleus. The source sets `9.4 T`, quadrupolar tensor eigenvalues `[-1e3 -2e3 3e3] Hz`, and Euler angles `[0 0 0]`. These are model inputs rather than experimental measurements; the tensor eigenvalues sum to zero. Spinach documents quadrupolar interaction tensors in Hz in its [g2spinach knowledge page](../../interfaces/g2spinach.md).

## Fokker–Planck label, actual call, and acquisition

The source header identifies Fokker–Planck theory and says perturbative corrections to the rotating-frame transformation are not applied. The code's actual simulation call is `singlerot(spin_system,@acquire,parameters,'nmr')`, not `floquet` or `gridfree`; the page therefore distinguishes the stated formalism from the invoked routine. The source estimates seconds for calculation time, not a measured runtime.

The rotor axis is `[1 1 1]` and rate `1000 Hz`; powder settings are `leb_2ang_rank_17` and `max_rank=17`. Acquisition uses a `2e4 Hz` sweep, 512 points, zero-fill 4096, offset 0, and ppm axis units; it selects `2H`, leaves `decouple={}`, and sets `invert_axis=1`. Both initial state and receiver are `L+` on `2H`.

The returned FID is exponentially apodised with parameter `6`, Fourier transformed with `fftshift(fft(fid,parameters.zerofill))`, and the real spectrum is plotted using `plot_1d`. This is a computed spectrum; the source does not provide experimental measured output or a numerical comparison.
