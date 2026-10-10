# kernel/integrity/patrol.m

Source: [MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/integrity/patrol.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=patrol.m).

- Signature: `patrol(test_subject)`

## Purpose

Checks MATLAB syntax and runs selected example files. The source comment describes its use as a continuous server-side safeguard, but one call to `patrol` performs a finite pass through the selected examples; the function itself does not contain a perpetual service loop.

## Selection and execution

The routine recursively lists `.m` files under `examples`. With the default empty subject it selects all listed files. A nonempty character-string subject selects a file when the text occurs either in the file's pathname or in one of its lines. The exception list is empty in the source.

It shuffles MATLAB's random-number generator, then repeatedly chooses a remaining selected file at random. For each file it checks `checkcode`; any diagnostic opens that file in the editor and raises an error. Otherwise it changes to the file's directory and evaluates the example, then flushes the display and pauses for one second before proceeding. The example's return value is not captured by `patrol`. The routine does not restore the prior working directory or random-number-generator state in this code.

## Inputs, outputs, and units

`test_subject` defaults to `''` when omitted and must be a character array; other supplied types fail in `grumble(test_subject)`. The function returns no value. Example execution may have its own effects and outputs, but `patrol` does not define a physical equation, normalisation, matrix output shape, or Hz/angular-frequency convention.

## Source guard

The selection is driven by literal `contains` checks against both file contents and the full pathname. Syntax diagnostics are fail-fast: the editor is opened for the affected file and the call stops with an error.

Related checks: [`existentials.m`](./existentials.md) checks startup prerequisites and path visibility; [`exorcise.m`](./exorcise.md) checks source conventions.
