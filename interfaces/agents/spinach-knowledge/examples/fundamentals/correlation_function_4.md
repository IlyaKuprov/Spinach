# examples/fundamentals/correlation_function_4.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/correlation_function_4.m`
- Signature: `correlation_function_4()`
- Total lines: 96

## Purpose

Computes rotational correlation functions using a Monte-Carlo method and compares them to the analytical results returned by Spinach kernel for the following correlation function: G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'> The sigma parameters refer to the rates of rotation and the four indices to the Wigner functions being correlated. High-rank isotropic rotational diffusion tested here. Calculation time: minutes.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Computes rotational correlation functions using a Monte-Carlo method and
- compares them to the analytical results returned by Spinach kernel for
- the following correlation function:
- G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'>
- The sigma parameters refer to the rates of rotation and the four indices
- to the Wigner functions being correlated. High-rank isotropic rotational
- diffusion tested here.
- Calculation time: minutes.
- Set testing parameters
- Convert indices from [-L,L] to [1,2*L+1]
- % Numerical Monte-Carlo calculation
- Number of points and lags

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `DCMT()`, `angles()`, `squeeze()`, `dcm2euler()`, `wigner()`, `xcorr()`, `ifftshift()`, `create()`, `basis()`, `corrfun()`, `kfigure()`, `lags()`, `cf_mc()`, `xlim()`, `kylabel()`, `kxlabel()`.
