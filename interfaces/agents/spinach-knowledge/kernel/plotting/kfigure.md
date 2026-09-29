# kernel/plotting/kfigure.m

Source: [kernel/plotting/kfigure.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/kfigure.m) · Source-listed Wiki URL: [Wiki](https://spindynamics.org/wiki/index.php?title=pauli.m)

- Signature: `handle = kfigure(varargin)`

## Behaviour and figure defaults

Before creating a figure, `kfigure` sets these four defaults on MATLAB's root object `groot`:

- `DefaultFigurePosition`: `[680 458 560 420]`
- `DefaultFigureWindowStyle`: `normal`
- `DefaultFigureMenuBar`: `figure`
- `DefaultFigureToolbar`: `figure`

The position value is retained as supplied to the root default; the function does not set a figure coordinate unit. The source describes these as the pre-R2025a settings. It then calls `figure(varargin{:})` and returns that handle. The root defaults are global MATLAB figure defaults and remain changed after the call, affecting later figures unless reset elsewhere.

All input arguments are forwarded to MATLAB's `figure`; this function performs no input validation. It creates a figure but does not calculate plot data, axes limits, or physical quantities.
