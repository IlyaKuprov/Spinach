# tests/lib/test_close.m

- Signature: `result=test_close(result,label,observed,reference,abs_tol,rel_tol,why)`

## Purpose

Adds a numerical regression check with tolerances and explanation. Syntax: result=test_close(result,label,observed,reference,abs_tol,rel_tol,why)

## Physical / mathematical content

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
