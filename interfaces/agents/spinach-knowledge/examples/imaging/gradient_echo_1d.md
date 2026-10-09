# examples/imaging/gradient_echo_1d.m

A one-dimensional gradient-echo simulation with diffusion and flow; the source labels its calculation time as seconds and credits Ahmed Allami and Ilya Kuprov. [Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/gradient_echo_1d.m)

## Model and sequence

The system contains one `1H` spin, `sys.magnet=5.9` and scalar chemical shift `1.0`. The source does not annotate units for those two values. It uses the `sphten-liouv` formalism with no basis approximation. The sequence calls `imaging(spin_system,@grad_echo,parameters)`: offset `0.0`, gradient amplitude `5e-6`, step duration `2e-4`, and 200 steps. No unit is attached to those gradient parameters in this example.

The one-dimensional geometry is `dims=0.30` with 100 points and `{'period',3}` differentiation. The relaxation phantom and operator are empty. Initial and coil spatial phantoms are uniform ones, with `Lz` for the initial state and `L+` for detection. Flow is set to `u=ones(100,1)` and diffusion to `1e-6`; the source does not state units for the geometry, flow, or diffusion values.

## Output and interpretation

The returned echo is plotted as its real part against a 401-point axis from `-0.04` to `+0.04` seconds (the range is calculated from the configured step duration and step count). The plot labels signal intensity in arbitrary units. The example does not define an RF-pulse amplitude or duration separately; avoid attributing either to this sequence from this file alone.
