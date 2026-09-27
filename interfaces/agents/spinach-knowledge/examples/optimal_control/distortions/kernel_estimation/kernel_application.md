# examples/optimal_control/distortions/kernel_estimation/kernel_application.m

- Signature: `kernel_application()`

## Purpose

HiPER instrument filter function kernel application to a complicated shaped pulse and a comparison with expe- rimental measurement at the instrument. Rob Hunter, Hassane el-Mkami, Graham Smith, Yujie Zhao, Shebha Anandhi Jegadeesan, Guinevere Mathies, Ilya Kuprov

## Physical / mathematical content

- The script combines the input pulse’s real and imaginary parts into a complex waveform, convolves it with the HiPER antenna kernel, and truncates the result to the input pulse length for comparison with the measured output.

## Numerical / algorithmic content

## Implementation structure

- HiPER instrument filter function kernel application to
- a complicated shaped pulse and a comparison with expe-
- rimental measurement at the instrument.
- Rob Hunter, Hassane el-Mkami, Graham Smith,
- Yujie Zhao, Shebha Anandhi Jegadeesan,
- Guinevere Mathies, Ilya Kuprov
- Read the input pulse
- Plot the input pulse
- Load appropriate HiPER kernel
- Compute and plot the convolution
- Read the measured pulse
- Plot the measured pulse
