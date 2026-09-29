# examples/fundamentals/perturb_theory.m

- Signature: `perturb_theory()`
- Source: [examples/fundamentals/perturb_theory.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/perturb_theory.m)

## Purpose

Compare Rayleigh-Schrodinger perturbation theory (RSPT) and Van Vleck perturbation theory (VVPT) for the energies and eigenvectors of the same finite Hermitian perturbation problem. The source describes the eigenvector representations as different while the energies agree; the example compares both perturbative results with direct diagonalisation.

## Model and assumptions

The no-argument script uses a 512-dimensional `pauli(512)` operator, with `H0=full(sigma.z)` as the source's Zeeman term. It generates a dense random complex matrix `X` and forms the Hermitian perturbation `H1=(1/25)*(X+X')/2`; there is no fixed random seed in the function. It is a numerical comparison on this constructed matrix, not a spin-system simulation or a statistical study.

## Checks and output

For each order 1 through 10, `rspert(diag(H0),H1,n)` and `vvpert(diag(H0),H1,n)` provide energies, which are sorted in descending order. The exact comparator is the sorted real spectrum of `H0+H1`. The plotted energy residual is the relative 2-norm of the difference between the exact energy vector and the zero-order or nth-order perturbative vector, divided by the exact energy-vector norm. A second plot compares eigenvector overlaps with the exact eigenvectors (sorted by their eigenvalues); RSPT supplies eigenvectors directly, whereas the VVPT output is a generator that the script exponentiates with `expm` before comparison. The overlap residual uses `norm(abs(V'*V_inf)-eye(size(V)),2)`, making it insensitive to eigenvector phase.

The script plots both residual sequences against orders 0 through 10 on logarithmic vertical axes. It does not specify a pass/fail tolerance, print a numerical table, or assert that either method meets a threshold; curve values depend on the random draw.

## Callable context

Run `perturb_theory()` with Spinach functions `pauli`, `rspert`, and `vvpert`, plus MATLAB linear algebra and plotting. The function takes no arguments and produces figures rather than returned values.
