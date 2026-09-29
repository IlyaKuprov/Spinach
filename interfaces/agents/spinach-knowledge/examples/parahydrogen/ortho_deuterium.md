# examples/parahydrogen/ortho_deuterium.m

- Signature: `ortho_deuterium()`
- Source: [examples/parahydrogen/ortho_deuterium.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/parahydrogen/ortho_deuterium.m)

## Purpose

Simulates the ortho-deuteration spectrum of acrylonitrile shown in Figure 1 of Natterer, Greve, and Bargon, cited by the source as [doi:10.1016/S0009-2614(98)00784-2](https://doi.org/10.1016/S0009-2614(98)00784-2). The source estimates a simulation time of seconds.

## Spin model and sequence

This is a five-spin product model: three `1H` spins and two `2H` spins, with the deuterons at indices 3 and 5. The code labels the setup “continuous deuteration,” sets `options.dephasing=1`, and calls `deut_pair` to obtain the deuteron-pair singlet `S` and quintet components `Q`. The starting operator is exactly `S+Q{1}+Q{2}+Q{3}+Q{4}+Q{5}`; the code does not select only the singlet or specify a parahydrogen singlet. There is no explicit reaction-rate, catalyst-exchange, or chemical-kinetics propagation in this script; it acquires a spectrum from the configured product spin system.

The field setting is `4.697` T, labelled “Bargon's magnet.” Under Spinach's nuclear Zeeman and scalar-coupling unit conventions, the chemical shifts are `{1.2, 1.2, 1.2, 2.3, 2.3}` ppm in isotope order `{1H, 1H, 2H, 1H, 2H}`. The source assigns scalar couplings `J(1,3)=2.0`, `J(2,3)=2.0`, `J(4,5)=2.0`, `J(3,4)=1.2`, `J(1,5)=1.2`, `J(2,5)=1.2`, and `J(3,5)=0.2` Hz; other matrix entries are not explicitly assigned. It uses the Zeeman Hilbert-space formalism with no basis approximation.

For the acquisition, the detected spins and coil are `2H`; the pulse operator is `Ly` on `2H` with angle `pi/4`. The settings are a 50 Hz transmitter offset, 120 Hz sweep, 1024 points, and 4096-point zero filling; the plotted axis is in ppm. The FID receives exponential apodisation with parameter 6, is Fourier transformed, and its real spectrum is plotted.

## Interpretation boundary

This is a simulated deuterium spectrum, not an experimental dataset. The source's `options.dephasing=1` is retained as a code setting; the script does not state that it is an exchange rate or a measured dephasing constant. Its placement in the parahydrogen examples directory does not make the modeled initial state a parahydrogen singlet.
