# kernel/integrity/patrol.m

- Signature: `patrol(test_subject)`

## Purpose

Runs a continuous patrol over Spinach example files to catch problems during development.

## Physical / mathematical content

This is a development and example-integrity utility; it does not define a physical model.

## Numerical / algorithmic content

The function repeatedly selects an example at random, checks it with Matlab's `checkcode`, and runs it unless it is listed as an exception. It waits one second between iterations and continues until interrupted.

## Parameters / inputs

- `test_subject` — a character string. When nonempty, only example files whose contents or path contain this string are included; an empty value selects all examples.

## Outputs

The output is whatever the individual examples display or return while they run. The patrol reports the number of selected files.

## Implementation structure

After filtering the example tree, the function errors if no files match. It shuffles the random-number generator, then loops over the selected files by choosing one at random, running the syntax check, changing to that file's directory, and evaluating the example.
