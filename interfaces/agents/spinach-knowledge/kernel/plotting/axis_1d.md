# kernel/plotting/axis_1d.m

Source: [kernel/plotting/axis_1d.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/axis_1d.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=axis_1d.m)

- Signature: `[ax,ax_label]=axis_1d(spin_system,parameters)`

## Purpose

Build the coordinate vector and text label for a one-dimensional spectrum axis. It returns axis data; it does not draw a plot, choose colours, or normalise a spectrum.

## Axis construction

`parameters.sweep` is either one real numeric value or two. For a scalar sweep, `parameters.offset` must be present and scalar; the function delegates construction to `ft_axis(parameters.offset,parameters.sweep,parameters.zerofill)`. The offset is relative to the magnet frequency in Hz. For a two-element sweep, the first endpoint must be less than the second, `offset` must be absent, and the Hz coordinates are `linspace(sweep(1),sweep(2),zerofill)`. Both paths produce the requested number of points and report the endpoint range through `report`.

The spin's reference frequency is calculated as `basefrq=-spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`. Any characters before the first letter in the spin label are formatted as an isotope superscript in the returned label when present; for example, the existing input example `{'1H'}` is accepted as the one-element `parameters.spins` cell.

## Coordinate units

The following conversions act on the Hz axis `ax`; `magnet` denotes `spin_system.inter.magnet`, and `freeg` denotes `spin_system.tols.freeg`.

| `parameters.axis_units` | Returned coordinates | Label meaning |
| --- | --- | --- |
| `'ppm'` | `-1e6*ax/basefrq` | isotope chemical shift in ppm |
| `'Gauss'` | `1e4*(magnet-2*pi*ax/spin('E'))` | magnetic induction in Gauss |
| `'mT'` | `1e3*(magnet-2*pi*ax/spin('E'))` | magnetic induction in mT |
| `'Hz'` | unchanged | isotope offset frequency in Hz |
| `'kHz'` | `1e-3*ax` | isotope offset frequency in kHz |
| `'MHz'` | `1e-6*ax` | isotope offset frequency in MHz |
| `'GHz'` | `1e-9*ax` | isotope offset frequency in GHz |
| `'MHz-labframe'` | `-1e-6*(ax-basefrq)` | isotope frequency in MHz |
| `'GHz-labframe'` | `-1e-9*(ax-basefrq)` | isotope frequency in GHz |
| `'gtensor'` | `-freeg*(ax-basefrq)/basefrq` | g-tensor coordinate scaled by the free-electron g value |
| `'points'` | `1:parameters.zerofill` | digitisation-point index |

For magnetic-field conversion the source uses the electron spin frequency through `spin('E')`. The label is returned as text for use by the caller; isotope superscripting is label formatting, not a graphics operation performed here.

## Inputs and outputs

- `parameters.sweep`: sweep width (scalar) or ascending Hz endpoints (two values).
- `parameters.offset`: required scalar only with a scalar sweep.
- `parameters.zerofill`: positive integer axis length.
- `parameters.axis_units`: one of the unit names in the table.
- `parameters.spins`: one-element cell naming the spin, such as `{'1H'}`.
- `ax`: row vector of axis coordinates; `ax_label`: axis-label text.
