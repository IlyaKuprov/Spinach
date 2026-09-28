# kernel/integrity/sniff.m

- Signature: `sniff(action)`

## Purpose

Compares Spinach `.m` files with the integrity baseline recorded in `smells.mat` by `rearm`. It reports files whose names or contents no longer match that baseline.

## Physical / mathematical content

This is a source-integrity utility; it does not model a physical system.

## Numerical / algorithmic content

For each included file, the routine hashes its filename and contents, then checks whether the resulting value is present in the saved `smells` list.

## Parameters / inputs

- `action` — `'none'` prints the names of flagged files; `'open'` opens them in the editor. The default is `'none'`.

## Outputs

No return value. Flagged files are reported as `smells fishy`; if all checks pass, the function prints the comment and code line counts and an all-clear message.

## Implementation structure

The function loads `smells.mat` and scans `.m` files under `kernel`, `interfaces`, `experiments`, and `etc`, applying the exception list. It ignores blank and comment lines when counting code and comment lines.
