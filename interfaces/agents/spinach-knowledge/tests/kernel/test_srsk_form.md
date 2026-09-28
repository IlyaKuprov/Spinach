# tests/kernel/test_srsk_form.m

- Signature: `result=test_srsk_form()`

## Purpose

Tests SRSK formalism restrictions and supported source relaxation.

## Physical / mathematical content

A rapidly relaxing 14N source broadens its scalar-coupled proton in spherical-tensor Liouville space. The test checks the Abragam longitudinal and transverse rates of the additive SRSK contribution.

## Numerical / algorithmic content

Checks that Zeeman Liouville SRSK is explicitly refused across lab-frame, secular, and diagonal retention and zero, IME, and dibari equilibrium choices, including when extended T1/T2 is requested. Separately checks a Zeeman Lindblad control without SRSK and the supported spherical-tensor Liouville SRSK contribution on longitudinal and complex transverse proton states.

## Syntax

`result=test_srsk_form()`

## Parameters / inputs

None. The test constructs its own spin-system fixtures.

## Outputs

`result` is a regression test result with explanatory messages.

## Header notes

Contact: talos@spindynamics.org.