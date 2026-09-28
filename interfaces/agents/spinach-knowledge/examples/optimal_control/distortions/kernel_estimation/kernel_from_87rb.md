# examples/optimal_control/distortions/kernel_estimation/kernel_from_87rb.m

- Signature: `kernel_from_87rb()`

## Purpose

Transmitter and probe distortion kernel of a 400 MHz Bruker spectrometer fitted with a 4 mm Phoenix MAS probe, estimated from an oscilloscope recording of an eight-block XiX waveform on the 87Rb channel. The wall clock record is heterodyned into the rotating frame, the carrier frequency is refined from the residual phase drift accumulated within each XiX half-block, the remaining constant phase is removed, and the FIR kernel is obtained by linear least squares from the ideal and measured waveforms.

## Physical / mathematical content

- The ideal XiX waveform alternates sign every half-block and is zero after the eighth block. A 40-tap FIR kernel is estimated from the ideal and measured waveforms and normalised to unit DC response.

## Numerical / algorithmic content

- The oscilloscope record is heterodyned at a refined carrier frequency, phase-corrected, amplitude-normalised, and resampled onto a 100 ns kernel grid. The kernel and its zero-filled Fourier magnitude response are plotted.

## Implementation structure

- Transmitter and probe distortion kernel of a 400 MHz Bruker
- spectrometer fitted with a 4 mm Phoenix MAS probe, estimated
- from an oscilloscope recording of an eight-block XiX wave-
- form on the 87Rb channel.
- The wall clock record is heterodyned into the rotating frame,
- the carrier frequency is refined from the residual phase drift
- accumulated within each XiX half-block, the remaining constant
- phase is removed, and the FIR kernel is then obtained by line-
- ar least squares from the ideal and the measured waveforms.
- The heterodyne is an analytic signal demodulation, which is
- zero-phase: the rotating frame signal is not delayed with re-
- spect to the oscilloscope record.
