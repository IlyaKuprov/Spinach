# interfaces/hiper/spinach2hiper.m

[Canonical source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/hiper/spinach2hiper.m) · [Wiki page](https://spindynamics.org/wiki/index.php?title=spinach2hiper.m)

**Call:** `spinach2hiper(file_name,amp,phi,off,dt)`. The function writes `[file_name '.csv']` and has no return value.

`file_name` is required to be a character string; the documented convention is a basename without an extension, and the implementation always appends `.csv`. `amp` and `phi` must be equal-length, finite, real numeric vectors. Amplitudes are written unchanged; phases are input in radians, converted to degrees, and wrapped to [0, 360]. `off` is a finite real numeric scalar in Hz and is repeated for each slice after conversion to MHz by multiplying by `1e-6`. `dt` is a positive finite real numeric scalar in seconds, converted to nanoseconds by multiplying by `1e9`.

The CSV table columns are `time_ns`, `freq_MHz`, `phase_deg`, and `amplitude`. Slice start times are `cumsum(dt_ns)-dt_ns`, so the first row starts at zero and subsequent starts are separated by one slice duration. The input vectors are columnised for the table. The source cautions that instrument phase direction may differ; test positive and negative phases as appropriate for the instrument.

The exporter calls `wrapTo360`, `table`, and `writetable`; it performs no amplitude scaling or hardware calibration.
