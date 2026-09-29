# examples/fundamentals/high_spin_system_1.m

- MATLAB implementation: [examples/fundamentals/high_spin_system_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/high_spin_system_1.m)

- Signature: `high_spin_system_1()`

## Physical question and spin model

This pulse-acquire NMR example illustrates the expected splitting of proton lines from a hypothetical scalar coupling to 235U. It builds a four-spin system at `14.1` T with isotopes `1H`, `235U`, `1H`, `1H`; the source assigns scalar shifts `-0.5`, `0.0`, `2.5`, and `1.3` ppm, respectively. The basis is `sphten-liouv` with approximation `none`. The specified scalar couplings are `J12=100` Hz and `J34=50` Hz, with `J44=0` also assigned.

## Acquisition and processing

The acquisition observes `1H`, starts from the `L+` state for `1H`, and uses the corresponding `L+` state as the coil. Decoupling is empty and the offset is 0. The sweep width is `3500` Hz, with `1024` acquired points and zero filling to `4096`; the displayed axis is in ppm and inverted. The code simulates an NMR FID with `liquid(spin_system,@acquire,parameters,'nmr')`, applies exponential apodisation with parameter 6, computes `fftshift(fft(fid,4096))`, and plots the real spectrum using `plot_1d`.

The source comment states the expected qualitative outcome, namely splitting from the hypothetical uranium coupling; the script does not encode a numeric peak-position check or report a measured spectrum. Exact plotted output depends on the simulation and processing settings above. Source: [examples/fundamentals/high_spin_system_1.m](../../../../../examples/fundamentals/high_spin_system_1.m).