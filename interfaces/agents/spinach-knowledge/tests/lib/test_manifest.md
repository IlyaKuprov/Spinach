# tests/lib/test_manifest.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/lib/test_manifest.m`
- Signature: `manifest=test_manifest()`
- Total lines: 136

## Purpose

Returns the registry of hand-written Spinach regression tests for the test runner.

## Physical / mathematical content

The manifest does not perform a physical calculation; its entries identify tests covering physical models, numerical methods, interfaces, and utility functions.

## Numerical / algorithmic content

Each entry stores a stable test ID, a descriptive name, and the function name to run. The cooperative-gradient, rotor-assumption, and tensor-train phase tests are registered. This catalogue does not execute the tests.

## Outputs

- `manifest` — structure array with stable `id`, descriptive `name`, and callable `function` fields.

