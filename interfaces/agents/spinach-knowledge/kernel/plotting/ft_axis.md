# kernel/plotting/ft_axis.m

- Signature: `ax=ft_axis(offset,sweep,npoints)`

## Purpose

Returns a row vector of Fourier-frequency coordinates with spacing `sweep/npoints`. It is an axis-construction helper, not a Fourier transform.

## Axis construction

The code first forms `linspace(-sweep/2,sweep/2,npoints+1)`. For odd `npoints` it drops the first value, shifts the remaining values left by half a bin, and adds `offset`; for even `npoints` it drops the last value and adds `offset`. Thus odd point counts give bins symmetric about the offset; even counts include the lower edge `offset-sweep/2` and stop one bin below the upper periodic edge `offset+sweep/2`. The coordinate units are those supplied for `offset` and `sweep`.

## Inputs and guards

- `offset` - real numeric scalar centre frequency.
- `sweep` - positive real numeric scalar frequency span.
- `npoints` - real numeric integer greater than 2.

Invalid inputs raise an error. The output `ax` is a 1-by-`npoints` row vector.

## Links

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/ft_axis.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ft_axis.m)
