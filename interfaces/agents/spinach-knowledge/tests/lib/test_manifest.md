# tests/lib/test_manifest.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/lib/test_manifest.m`
- Signature: `manifest=test_manifest()`
- Total lines: 134

## Purpose

Registers Spinach's hand-written regression tests for the test runner.

## Physical / mathematical content

The registry covers physical models, numerical methods, interfaces, and utility functions; the calculations reside in the registered tests.

## Numerical / algorithmic content

Each entry associates a stable test identifier and a descriptive name with a test function name. The merged registry includes the branch-specific regression and the cooperative-gradient and rotor-assumption tests from main; this function does not execute them.

## Syntax

```matlab
manifest=test_manifest()
```

## Parameters / inputs

None.

## Outputs

`manifest` is a structure array with `id`, `name`, and `function` fields.

## Header notes

The manifest stores test metadata only.
