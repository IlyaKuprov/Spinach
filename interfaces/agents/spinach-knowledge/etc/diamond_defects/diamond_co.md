# etc/diamond_defects/diamond_co.m

- MATLAB implementation: [etc/diamond_defects/diamond_co.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_co.m)

- Signature: `[sys,inter]=diamond_co(parameters)`
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=diamond_co.m)
- Magnetic parameters: Nadolinny et al., *Crystals* **7**, 237 (2017), [https://doi.org/10.3390/cryst7080237](https://doi.org/10.3390/cryst7080237)
- Source author: alexey.bogdanov@weizmann.ac.il

## Purpose

Construct Spinach system and interaction structures for an electron–cobalt defect in diamond using the tabulated parameters for the `o4` or `nlo2` centre.

## Call and inputs

`[sys,inter]=diamond_co(parameters)` takes exactly one structure with both fields:

- `parameters.centre` — character string `'o4'` or `'nlo2'`; selection is case-insensitive.
- `parameters.orientation` — character string `'111'`, `'110'`, or `'100'`, selecting the crystal plane normal aligned with the magnetic field. The code matches these values exactly.

## Magnetic tensors and units

The centre selects principal values:

| Centre | Electron `g` principal values | `59Co` hyperfine principal values (mT) | Frame angle `α` |
|---|---|---|---:|
| `o4` | [2.3463, 1.8438, 1.7045] | [8.86, 6.43, 5.82] | 29° |
| `nlo2` | [2.3277, 1.7982, 1.7149] | [8.24, 6.57, 5.76] | 28° |

The source constructs a principal-axis frame from `y=[0,1,1]/√2`, rotates the x axis by `α` (degrees), orthogonalises the frame, and enforces right-handedness. It builds the electron tensor as `frame·diag(g)·frameᵀ`; the cobalt hyperfine tensor is built analogously after converting the tabulated mT values to Spinach frequency units using `abs(spin('E'))/(2π)·10⁻³`. A rotation aligns the selected crystallographic normal with the magnetic-field z axis, and is applied to both tensors.

## Outputs

- `sys` — isotopes `{'E','59Co'}`.
- `inter` — electron Zeeman tensor at `inter.zeeman.matrix{1}` and electron–cobalt hyperfine tensor at `inter.coupling.matrix{1,2}`.

The routine specifies the electron Zeeman and electron–nuclear hyperfine tensors; it does not add a cobalt Zeeman tensor here. Unsupported centre or orientation strings raise an error.
