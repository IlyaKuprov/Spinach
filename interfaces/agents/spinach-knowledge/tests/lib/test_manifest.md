# tests/lib/test_manifest.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/lib/test_manifest.m`
- Signature: `manifest=test_manifest()`
- Total lines: 132

## Purpose

Returns the metadata that registers Spinach's hand-written, physically motivated regression tests.

## Physical / mathematical content

The registry spans algebra, propagation, spin-system construction, spectroscopy, and numerical utilities; the scientific calculations reside in the registered test functions, not in this metadata function.

## Numerical / algorithmic content

Each entry associates a stable test identifier and a descriptive name with a test function name. `run_tests` uses this list to select and execute tests. The stack-reduction regression is registered as `kernel/stack_reduce`, implemented by `test_stack_reduce`.

## Syntax

```matlab
manifest=test_manifest()
```

## Parameters / inputs

None.

## Outputs

`manifest` is a structure array with fields `id`, `name`, and `function`, containing test identifiers, descriptions, and function names, respectively.

## Header notes

The manifest stores metadata only; it does not execute the registered tests.
