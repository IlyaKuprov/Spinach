# tests/lib/test_true.m

- Signature: `result=test_true(result,label,condition,why)`

## Purpose

Adds a logical regression check with a clear message. Syntax: result=test_true(result,label,condition,why)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- result -test result structure
- label -check label
- condition -logical pass/fail condition
- why -explanation of the right answer

## Outputs

- result -updated test result structure

## Implementation structure

- Adds a logical regression check with a clear message. Syntax:
- result=test_true(result,label,condition,why)
- result -test result structure
- label -check label
- condition -logical pass/fail condition
- why -explanation of the right answer
- result -updated test result structure
- Check the condition
- Record the pass message
