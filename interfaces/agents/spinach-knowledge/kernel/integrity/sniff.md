# kernel/integrity/sniff.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/integrity/sniff.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=sniff.m)

`sniff(action)` compares current Spinach .m files with the fingerprint list made by [`rearm`](rearm.md). It is an integrity check, not a model of a physical system.

The scan covers .m files recursively below `kernel`, `interfaces`, `experiments`, and `etc`, using the same basename-plus-content fingerprint as `rearm`. A fingerprint absent from `smells` is printed as `smells fishy: <path>`; with action `'open'`, the file is also opened in the editor. If all scanned fingerprints are found, `sniff` prints comment-line and code-line counts and an all-clear message. Blank lines are excluded from those counts, and only lines whose first character is `%` are counted as comment lines. These counting filters are applied after fingerprinting, so comments and blank lines still affect the fingerprint.

The saved values are tested by membership, not matched to paths. A changed file can be flagged, but an identical basename and content at another location can match; deleted files are not detected by scanning. `sniff` loads `smells.mat` by that bare filename and does not itself rebuild the baseline; call [`rearm`](rearm.md) to establish a new one.

## Inputs and outputs

`action` is an optional character input: omitted means `'none'`; the only accepted values are `'none'` and `'open'`. Other values raise an error. The function has no return value. Failure to load `smells.mat` is not caught by a source-level recovery guard.
