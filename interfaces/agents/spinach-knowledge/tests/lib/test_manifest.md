# tests/lib/test_manifest.m

- Signature: `manifest=test_manifest()`

## Purpose

Registers the regression tests available to the Spinach test runner. Adding a test file does not register it: its stable identifier, descriptive name, and callable function name must also appear in this manifest.

## Physical / mathematical content

None is evaluated here. The listed tests cover physical and numerical subjects, but the manifest itself is metadata, not an implementation of those theories or algorithms.

## Numerical / algorithmic content

Returns an ordered structure array with the same three fields for every registration. It does not run tests, calculate tolerances, inspect results, or infer registrations from files on disk. The runner uses these records to identify and invoke test functions.

## Syntax

`manifest=test_manifest()`

## Parameters / inputs

None.

## Outputs

`manifest` is a structure array. Each element contains `id` (stable test identifier), `name` (human-readable description), and `function` (test function name as a character vector).

## Header notes

The function returns regression-test metadata only.
