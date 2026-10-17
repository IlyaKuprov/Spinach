# kernel/plotting/plot_uf.m

- Signature: `plot_uf(spin_system,spectrum_uf,parameters)`

## Purpose

Plot an ultrafast constant-time 2D spectrum with a UF axis and a conventional axis. The routine constructs plotting coordinates and draws contours; it does not calculate the spectrum.

## Axis construction

Let `Ta=parameters.deltat*parameters.npoints` seconds. The conventional sweep is `sweep_conv=1/(2*Ta)` Hz, and the conventional axis is `sweep_conv*(-floor(nloops/2):(ceil(nloops/2)-1))/nloops+offset(2)` Hz. The first entry in `parameters.spins` supplies `gamma=spin(spins{1})` in rad/s/T; the code sets `k_max=gamma*Ga*Ta/(2*pi)` in m⁻¹, `t_max=2*Te`, and `constant_c=(-2*t_max)/dims`. It then uses `res_uf=abs(1/(constant_c*dims))` Hz and `sweep_uf=abs(k_max/constant_c)` Hz. The UF-axis sample count is `round(dims*k_max)`, which must equal `size(spectrum_uf,1)`; its values are `-sweep_uf/2+res_uf*(0:(uf_dim_size-1))+offset(1)`.

For `axis_units='ppm'`, the two axes are converted with `1e6*(2*pi*frequency)/(spin(spins{i})*spin_system.inter.magnet)`; for `'Hz'`, the frequency values are unchanged. The code then adds `offset_uf_cov` to the UF/F1 axis in the selected units. Offsets `offset(1)` and `offset(2)` are in Hz before conversion; `offset_uf_cov` is documented in the selected axis units. The constant-time-axis construction cites *Progress in Nuclear Magnetic Resonance Spectroscopy* **57** (2010), 241; no DOI is given in the source.

## Rendering and checks

The contour call uses `axis_f2` for X, `axis_f1` for Y, and `flipud(spectrum_uf)` for the matrix. It boxes and grids the current axes, labels X as `1Q / ppm` and Y as `MQ / ppm`, and reverses both axis directions. Those labels are literal source strings; the code does not compose them from `axis_units`. No explicit axis limits or colormap are set here, and the function returns no value.

The input must be a real numeric matrix. The parameter structure must supply spins, two offsets, `offset_uf_cov`, sample dimension, acquisition time step and point count, loop count, gradient amplitude, echo time, and axis units; the source checks these scalar values for the required finite/positive forms and accepts only `'ppm'` or `'Hz'`. A single spin entry is duplicated for both dimensions. The UF row count is checked against the rounded `dims*k_max` value.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/plot_uf.m)
- [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=plot_uf.m)
