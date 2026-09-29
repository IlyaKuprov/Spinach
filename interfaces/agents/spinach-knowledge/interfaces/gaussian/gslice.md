# interfaces/gaussian/gslice.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/gaussian/gslice.m) · [Spinach Wiki: gslice.m](https://spindynamics.org/wiki/index.php?title=gslice.m)

- Signature: `gslice()`

## Inputs and selection

This is an interactive MATLAB utility with no input arguments. It prompts for a Gaussian relaxed-geometry scan log and for a separate text file containing the Gaussian calculation header. It scans the log for the exact marker `-- Stationary point found.`; for each such marker it searches backwards for the nearest preceding `Standard orientation:` block.

From each selected orientation it reads the centre/atomic-number/coordinate rows and converts atomic numbers to element symbols. It writes a new input geometry for each selected scan minimum, using the provided header followed by the selected coordinates formatted to eight decimal places. The function does not compute or validate energies or property values at those geometries; it prepares input files only.

## Files and return behaviour

In the current working directory, it writes numbered `g16_input_N.gjf` files and a `compute.bat` batch file. The batch file invokes Gaussian 16 and sets `GAUSS_EXEDIR=C:\G16W` and `GAUSS_SCRDIR=C:\Temp`; both paths are hard-coded for Windows and must match the local installation. The function takes no arguments, returns no MATLAB output, and does not run the generated batch file. It relies on MATLAB's `uigetfile` dialogs and text/file I/O, plus its local periodic-table mapping.

[Spinach Wiki: gslice.m](https://spindynamics.org/wiki/index.php?title=gslice.m)
