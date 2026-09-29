# kernel/derivatives/fftdiff.m

Direct source: [kernel/derivatives/fftdiff.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/derivatives/fftdiff.m)
Spin Dynamics Wiki: [fftdiff.m](https://spindynamics.org/wiki/index.php?title=fftdiff.m)

## Purpose and interface

kern=fftdiff(order,npoints,dx) returns a length-npoints Fourier-domain multiplier for differentiating a uniformly sampled periodic real signal. Apply it to a signal with npoints entries in the Fourier-bin layout of fft:

    signal=[0 1 0 -1 0 1 0 -1];
    kern=fftdiff(1,8,0.25);
    derivative=real(ifft(fft(signal).*kern));

The kernel is a vector of spectral multipliers, not a matrix or a signal. Under the discrete Fourier representation, periodic Fourier modes are eigenvectors of differentiation; the corresponding multiplier for a mode with wave number k is (2*pi*1i*k/(npoints*dx))^order. The result is mapped back to sample space by the inverse FFT, and real keeps the real-valued result for real input signals.

## Frequency grid and units

The source uses centred integer mode indices and then ifftshift to arrange them in the bin order expected by fft. For odd npoints, the centred indices run from (1-npoints)/2 through (npoints-1)/2. For even npoints, they run from -npoints/2 through (npoints-1)/2, including the negative Nyquist bin. Multiplication by 2*pi/(npoints*dx) sets the angular wave-number scale; raising the imaginary multiplier to order gives the requested derivative order. The units are signal units divided by dx^order.

The method assumes periodic boundary conditions and uniform spacing dx. It is the Fourier/spectral alternative to the finite-stencil construction in [fdmat.m](fdmat.md).

## Input guards

The source requires order and npoints to be real numeric scalars, positive integers. It requires dx to be a real numeric scalar greater than zero. The error strings describe order and npoints as “non-negative,” but the actual predicates reject values below 1, so zero is not accepted.
