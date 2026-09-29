# kernel/plotting/ksgtitle.m

- MATLAB implementation: [kernel/plotting/ksgtitle.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/ksgtitle.m)

Apply the Spinach figure style to the overall title for the current plot grid.

## Call

`ksgtitle(x)`

Pass a MATLAB character vector for `x`. The function checks `ischar(x)`; other types, including MATLAB string scalars, raise `x must be a character string`. The text is inserted directly inside a LaTeX `\textbf{...}` wrapper, so it is rendered bold by the LaTeX interpreter and is not escaped or validated by this helper. It does not accept additional `sgtitle` options.

## Effect

The helper calls `sgtitle` with the wrapped text and `'Interpreter','latex'`, creating or updating the overall title for the current grid. It returns no output; the effect is on the graphics title.

There is no physical or numerical calculation here: this is a plotting-style convenience wrapper.

## Source link

[`ksgtitle.m` on the Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=ksgtitle.m)
