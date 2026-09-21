# tests/lib/test_close.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/lib/test_close.m`
- Signature: `result=test_close(result,label,observed,reference,abs_tol,rel_tol,why)`
- Total lines: 71

## Purpose

Adds a numerical regression check with tolerances and explanation. Syntax: result=test_close(result,label,observed,reference,abs_tol,rel_tol,why)

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

## Parameters / inputs

- result -test result structure
- label -check label
- observed -value produced by Spinach
- reference -independently known right answer
- abs_tol -absolute tolerance
- rel_tol -relative tolerance
- why -explanation of the right answer

## Outputs

- result -updated test result structure

## Implementation structure

- Adds a numerical regression check with tolerances and explanation. Syntax:
- result=test_close(result,label,observed,reference,abs_tol,rel_tol,why)
- result -test result structure
- label -check label
- observed -value produced by Spinach
- reference -independently known right answer
- abs_tol -absolute tolerance
- rel_tol -relative tolerance
- why -explanation of the right answer
- result -updated test result structure
- Convert sparse arrays for norm evaluation
- Check dimensions first

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isequal()`, `double()`, `observed()`, `reference()`, `any()`, `num2str()`.
