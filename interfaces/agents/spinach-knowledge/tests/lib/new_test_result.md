# tests/lib/new_test_result.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/lib/new_test_result.m`
- Signature: `result=new_test_result(id,name,purpose)`
- Total lines: 30

## Purpose

Creates a regression test result structure. Syntax: result=new_test_result(id,name,purpose)

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

## Parameters / inputs

- id -stable test identifier
- name -short human-readable test name
- purpose -one-sentence purpose statement

## Outputs

- result -test result structure

## Implementation structure

- Creates a regression test result structure. Syntax:
- result=new_test_result(id,name,purpose)
- id -stable test identifier
- name -short human-readable test name
- purpose -one-sentence purpose statement
- result -test result structure
- Build the result structure
