# kernel/plotting/ktitle.m

- MATLAB implementation: [kernel/plotting/ktitle.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/ktitle.m)

Set the title of the current MATLAB axes in the Spinach plotting style.

## Call

`ktitle(x)`

Pass a MATLAB character vector for `x`. The function checks `ischar(x)`; a non-character input raises `x must be a character string`. It inserts the text directly inside a LaTeX `\textbf{...}` wrapper and passes it to `title` with the LaTeX interpreter, so the title is bold. The text is not escaped or checked for valid LaTeX, and the wrapper offers no extra `title` property arguments.

## Effect

The call creates or updates the title of the current axes; it is not a figure/grid-wide title. It returns no output and performs no physical or numerical calculation.

## Source link

[`ktitle.m` on the Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=ktitle.m)
