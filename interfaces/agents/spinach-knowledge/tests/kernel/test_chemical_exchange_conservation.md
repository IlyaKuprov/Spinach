# tests/kernel/test_chemical_exchange_conservation.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_chemical_exchange_conservation.m`
- Signature: `result=test_chemical_exchange_conservation()`
- Total lines: 44

## Purpose

Tests conservation in two-site chemical exchange. Syntax: result=test_chemical_exchange_conservation()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test builds a symmetric two-site exchange model and checks that the
- kinetics generator conserves the total population over the two sites.

## Implementation structure

- Tests conservation in two-site chemical exchange. Syntax:
- result=test_chemical_exchange_conservation()
- result -regression test result with explanatory messages
- The test builds a symmetric two-site exchange model and checks that the
- kinetics generator conserves the total population over the two sites.
- Announce the test target
- State the kinetics target of the test
- Build a symmetric two-site exchange system
- Build the kinetics generator
- Closed Markov kinetics conserve total population by zero column sums

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `kinetics()`, `test_close()`.
