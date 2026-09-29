# kernel/utilities/apodisation.m

## Purpose

Apodises free induction decay (FID) data by applying window functions along each dimension, with optional halving of first points to satisfy Fourier transform symmetry requirements. Supports FIDs of any dimensionality.

Source: [apodisation.m on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/apodisation.m)

## Behaviour

- Validates inputs via an internal `grumble` subroutine: `fid` must be numeric, `winfuns` must be a cell array with one element per non-singleton dimension of `fid`, each element must be a cell array, and the window type must be one of the supported strings. Parameterised windows (`exp`, `gauss`, `kaiser`, `bad-z1`, `bad-z2`) require a finite real scalar parameter; all others take no parameters.
- If `fp_half` is not supplied, it defaults to `true`.
- Non-singleton dimensions are identified with `find(size(fid)>1)`; dimensions whose `winfuns` entry is empty are excluded as inactive.
- When `fp_half` is true, the first points along each active dimension are divided by 2 (indexing the hyperplane where that dimension equals 1), and a report is printed per dimension.
- For each active dimension, a window function vector of length `npts = size(fid,dim)` is built and applied by elementwise multiplication after reshaping to match the array dimensions.
- Supported window types:
  - `none`: all-ones window (no apodisation, but first-point halving still applies when enabled).
  - `crisp`: `cos(x).^8` half-bell with `x` from 0 to `pi/2`.
  - `exp`: `exp(-k*x)` with `x` from 0 to 1.
  - `gauss`: `exp(-k*(x.^2))` with `x` from 0 to 1.
  - `cos`: `cos(x)` half-bell with `x` from 0 to `pi/2`.
  - `sin`: `sin(x)` full bell with `x` from 0 to `pi`.
  - `sqcos`: `cos(x).^2` half-bell with `x` from 0 to `pi/2`.
  - `sqsin`: `sin(x).^2` full bell with `x` from 0 to `pi`.
  - `kaiser`: MATLAB `kaiser(npts,k)` window; the peak is in the middle of the FID.
  - `bad-z1`: `sinc(x*k)` over `x` from 0 to 1 (excluding the first point, which is set to 1 to avoid the singularity); emulates a misset Z1 shim.
  - `bad-z2`: `(fresnelc(x)+1i*fresnels(x))/x` with `x = sqrt(t*k)` for `t` from 0 to 1 (first point set to 1); computed in a `parfor` loop; emulates a misset Z2 shim.
  - The source header gives illustrative dimensionless `k` values for 1H NMR at 600 MHz: `bad-z1` uses 10 and `bad-z2` uses 40; these are not defaults enforced by the function.
- An unsupported window type raises an error.
- After each window is applied, a report line is printed naming the dimension and window type.

## Inputs and outputs

**Syntax:** `fid = apodisation(spin_system, fid, winfuns, fp_half)`

**Inputs:**

- `spin_system` — spin system object, used for reporting.
- `fid` — the FID. A column vector for 1D data; for higher dimensions, a matrix (2D, 3D, etc.) with the time origin at the `(1,1)`, `(1,1,1)`, ... corner.
- `winfuns` — cell array of window function specifications, one per non-singleton dimension, in the format `{{spec},{spec},...}`. Each specification is one of the window types listed above; an empty element marks the dimension as inactive.
- `fp_half` — optional logical; set to `false` to disable first-point halving, needed when multiple window functions are applied to the same dimension sequentially.

**Outputs:**

- `fid` — the apodised free induction decay.

## References

- [apodisation.m — Spinach Wiki](https://spindynamics.org/wiki/index.php?title=apodisation.m)
- [apodisation.m source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/apodisation.m)
