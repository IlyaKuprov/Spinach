# examples/fundamentals/correlation_function_5.m

- Signature: `correlation_function_5()`

## Purpose

Computes the following rotational correlation function G(k,m,p,q)=<R(k,m)*R(p,q)> where R is the 3D Cartesian rotation matrix, using the Monte-Carlo method. Calculation time: minutes.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Computes the following rotational correlation function
- G(k,m,p,q)=<R(k,m)*R(p,q)>
- where R is the 3D Cartesian rotation matrix, using the
- Monte-Carlo method.
- Calculation time: minutes.
- Set testing parameters
- Set number of points
- Generate angle track
- Preallocate rotation matrix array
- Loop over Monte-Carlo steps
- Generate a random rotation
- Get Monte-Carlo correlation function
