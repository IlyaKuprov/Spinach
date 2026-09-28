# kernel/plotting/int_2d.m

- Signature: `int_2d(spin_system,spectrum,parameters,ncont,...`

## Purpose

Plots positive and/or negative contours of a two-dimensional spectrum with nonlinear level spacing, then integrates regions selected interactively or loaded from an interval file. Syntax: `int_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs,filename)`.

## Physical / mathematical content

## Numerical / algorithmic content

- Positive levels are `(delta(2)-delta(1))*smax*linspace(0,1,ncont).^k + smax*delta(1)`; negative levels are `(delta(4)-delta(3))*smin*linspace(0,1,ncont).^k + smin*delta(3)`, where `smax` and `smin` are taken from the spectrum. `k=1` gives linear spacing; `k>1` increases contour density near the baseline.
- If the interval file does not exist, the user selects integration regions and they are saved. If it exists, its regions are loaded and integrated using spline interpolation and `integral2` (relative and absolute tolerances `1e-3`).

## Parameters / inputs

- `spin_system` - spin system used for plotting
- `spectrum` - real matrix containing the 2D NMR spectrum
- `parameters.sweep` - one or two sweep widths, in Hz
- `parameters.spins` - cell array with one or two character strings specifying the working spins
- `parameters.offset` - one or two transmitter offsets, in Hz
- `parameters.axis_units` - axis units: `ppm`, `Hz`, or `Gauss`
- `ncont` - number of contours; 20 is a reasonable value
- `delta` - four fractions specifying lower and upper levels for positive and negative contours; starting value `[0.02 0.2 0.02 0.2]`
- `k` - curvature of contour spacing; 1 gives linear spacing, while values greater than 1 increase sampling density near the baseline; 2 is a reasonable value
- `ncol` - number of colours in the colour map; around 256 is suitable
- `m` - colour-map curvature; 1 gives a linear red/blue ramp for positive/negative contours; 6 is suitable for high contrast
- `signs` - `positive`, `negative`, or `both`
- `filename` - interval-file path; if it does not exist, selected intervals are saved there; if it exists, its intervals are loaded and integrated

## Outputs

- Creates a figure and reports integrals over the selected or loaded intervals; creates the interval file when needed.
