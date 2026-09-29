# kernel/utilities/overwound.m

**Source:** <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/overwound.m>

## Purpose

Checks if a Fokker-Planck state vector has any spatial frequencies that its spatial grid is dangerously close to misrepresenting due to insufficient point count.

## Behaviour

- Syntax: `overwound(rho,spc_dim,spn_dim)`.
- Calls the internal consistency-checking subfunction `grumble(rho,spc_dim,spn_dim)` before proceeding.
- Determines the stack size as `size(rho,2)` and reshapes `rho` into `[spn_dim spc_dim stack_size]`.
- Permutes the spin dimension out of the way with `permute(rho,[2 3 4 1 5])`.
- For each spatial dimension X, Y, Z that has more than one point:
  - Errors with `you must increase point count in X dimension.` (or Y/Z) if the point count in that dimension is fewer than 10.
  - Computes the FFT along that dimension with `fftshift(fft(...))`.
  - Permutes and reshapes the transformed array to pull the frequency dimension forward.
  - Plots the summed absolute spectrum against a frequency axis `linspace(-1,1,...)`, labelled `X spatial freq. relative to Nyquist limit` (or Y/Z), with y-label `population density, a.u.`, using `kfigure`, `kxlabel`, `kylabel`, and `kgrid`.
- The `grumble` subfunction enforces:
  - `rho` must be a numeric array.
  - `spc_dim` must be a vector of three positive integers.
  - `spn_dim` must be a positive integer scalar.
  - The number of rows in `rho` must equal `prod([spn_dim spc_dim])`.
- The header note advises that when running spatial dynamics such as diffusion and flow with finite difference derivative operators, the spatial grid point count should be set to several times the minimum Nyquist value.

## Inputs and outputs

**Inputs:**

- `rho` — Fokker-Planck state vector or a bookshelf stack thereof.
- `spc_dim` — spatial dimensions of the Fokker-Planck problem, `[X Y Z]`.
- `spn_dim` — spin dimension of the Fokker-Planck problem.

**Outputs:**

- Figures and diagnostic messages to the console.

## References

- Spinach wiki: <https://spindynamics.org/wiki/index.php?title=overwound.m>
