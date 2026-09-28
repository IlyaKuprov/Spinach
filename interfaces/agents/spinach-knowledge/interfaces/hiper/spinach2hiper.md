# interfaces/hiper/spinach2hiper.m

- Signature: `spinach2hiper(file_name,amp,phi,off,dt)`

Exports phase-modulated optimal-control waveforms for Graham Smith's HiPER instrument to `file_name.csv`.

`file_name` is a character string without the extension. `amp` and `phi` are equal-length finite real vectors (amplitudes and phases in radians); `off` is a finite real transmitter offset in Hz; and `dt` is a positive finite real slice duration in seconds.

The CSV columns are `time_ns` (slice start times from zero, spaced by `dt*10^9` ns), `freq_MHz` (`off*10^-6`), `phase_deg` (phase converted to degrees and wrapped to [0,360]), and `amplitude`.

[Spinach documentation](https://spindynamics.org/wiki/index.php?title=spinach2hiper.m)
