# tests/lib/new_test_result.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/lib/new_test_result.m`
- Signature: `result=new_test_result(id,name,purpose)`
- Total lines: 30

## Purpose

Creates a regression test result structure. Syntax: result=new_test_result(id,name,purpose)

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Build the result structure; implemented by `result.id=id`.

### Key state/data transformations

- Lines 22: computes `result.id` using `result.id=id`.
- Lines 23: computes `result.name` using `result.name=name`.
- Lines 24: computes `result.purpose` using `result.purpose=purpose`.
- Lines 25: computes `result.status` using `result.status='RUNNING'`.
- Lines 26: computes `result.elapsed` using `result.elapsed=0`.
- Lines 27: computes `result.messages` using `result.messages={}`.
- Lines 28: computes `result.error` using `result.error=''`.

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
