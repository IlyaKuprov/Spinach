# examples/nmr_paramag/gau_density.m

- Signature: `gau_density()`
- Source: [examples/nmr_paramag/gau_density.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_paramag/gau_density.m)

## Calculation

This is a forward, three-dimensional PCS-field calculation for an electron-density array, not a fit. It makes a 100-by-100-by-100 mesh on each axis from -10 to 10, sets `sigma=0.5`, and evaluates `probden = exp(-(X.^2+Y.^2+Z.^2)/(2*sigma))/sqrt((2*pi)^3*sigma^3)`. Those grid and density parameters have no units stated in the source.

The susceptibility tensor is `R*diag([-0.1,-0.2,0.3])*R'`, with `R=euler2dcm(pi/4,pi/5,pi/6)`; the tensor units are not stated. `kpcs(probden,chi,[-10 10 -10 10 -10 10],[0 0 0],'fft')` returns the PCS volume, which is displayed with `volplot`.

## Scope and omissions

No nuclear coordinates, experimental shifts, fitted parameters, spectrum, magnetic field, or temperature are specified. The source is a generic Gaussian-density example, not a carbonic-anhydrase case; it does not identify a protein residue or metal site. It displays a field volume rather than reporting a fixed numeric result.
