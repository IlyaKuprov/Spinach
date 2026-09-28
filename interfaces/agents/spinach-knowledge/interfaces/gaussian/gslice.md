# interfaces/gaussian/gslice.m

- Signature: `gslice()`

Prompts for a Gaussian relaxed-geometry scan log and a text file containing the Gaussian input header. For each `-- Stationary point found.` flag, it selects the nearest preceding `Standard orientation:` block and uses its atomic numbers and coordinates.

The function writes numbered `g16_input_N.gjf` files and a `compute.bat` script in the current directory. The batch file runs the inputs with Gaussian 16 and contains hard-coded Windows paths (`C:\G16W` and `C:\Temp`); edit these for the local installation. The function takes no arguments and returns no MATLAB outputs.

[Source documentation](https://spindynamics.org/wiki/index.php?title=gslice.m)
