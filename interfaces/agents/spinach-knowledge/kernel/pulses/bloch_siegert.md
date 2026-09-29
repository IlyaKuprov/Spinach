# kernel/pulses/bloch_siegert.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/bloch_siegert.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=bloch_siegert.m)

- Signature: `[ctrl_opers,ctrl_coefs]=bloch_siegert(spin_system,ctrl_opers,ctrl_coefs)`

## Inputs and checks

The two control arguments are cell arrays with one operator and one real coefficient vector per control channel. Coefficients are in rad/s; all vectors must have the same length. The spin-system control configuration must have Bloch-Siegert corrections enabled through `optimcon()` and provide the channel isotope, channel index, and carrier-frequency settings used to build the response operators.

## Augmentation and ordering

For each input channel, the function builds its Bloch-Siegert response operator from the configured isotope and carrier frequency. It appends the response operators after the original operators, preserving input-channel order. The matching appended coefficient vector is the elementwise square of that channel's original real amplitude vector. Thus the number of samples and their positions in each vector are retained; this routine does not add a time grid, change pulse duration, or apply a separate phase parameter. It leaves any quadrature/channel convention to the caller's input ordering.

The returned arrays are prepared for downstream `shaped_pulse_xy` / GRAPE use. This function returns augmented arrays; it does not write a waveform or export file, and the source specifies no file format.
