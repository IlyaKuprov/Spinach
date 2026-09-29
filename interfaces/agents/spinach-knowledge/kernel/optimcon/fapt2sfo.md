# kernel/optimcon/fapt2sfo.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/fapt2sfo.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=fapt2sfo.m)

- Signature: `[wave,dt,time_grid]=fapt2sfo(fapt,time_grid)`

## Purpose and units

Converts frequency-amplitude-phase-time pulse events to a two-row single-frequency-origin waveform for GRAPE. Each `fapt` cell contains a real five-element vector `[frequency, amplitude, phase, start_time, end_time]`: Hz, rad/s, radians, seconds, seconds. Amplitude must be nonnegative and end time must be strictly greater than start time. Events are active at grid ticks satisfying `start_time <= t <= end_time`; overlapping events add.

For an event at frequency `f`, amplitude `a`, and phase `phi`, the X and Y rows accumulate `a*cos(2*pi*f*t+phi)` and `a*sin(2*pi*f*t+phi)` respectively. This is the documented anticlockwise rotation convention. With a drift offset `2*pi*f*Lz`, an event at `f` is on resonance; reversing the Y sign gives the opposite sense and an offset of `2*f` for nonzero `f`. At `f=0`, the sign change reflects the nutation axis to `-phi`.

## Time grid and outputs

If `time_grid` is omitted, the function spans from 0 to the latest event end time. It chooses `dt_hN = 1/(4*max(abs(frequency)))`, sets `npts = ceil(end_time/dt_hN)+1`, and uses `linspace` to make the grid; the returned `dt` is the second grid tick. If every event frequency is zero, an explicit grid is required. If supplied, `time_grid` is checked to be a real numeric row vector and `dt` is empty.

- `wave`: `2 x numel(time_grid)` array, with X then Y components, in rad/s.
- `dt`: generated-grid step in seconds, or empty when the grid was supplied.
- `time_grid`: row vector of sample times in seconds.
