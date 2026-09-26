# tests/kernel/test_dynamic_parse_text_reporting_suite.m

- Signature: `result=test_dynamic_parse_text_reporting_suite()`

## Purpose

Tests deterministic parsing, text, and safe reporting utilities. Syntax: result=test_dynamic_parse_text_reporting_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks operator-specification parsing, isotope predicates,
- label lookup, silent reporting calls, and polyadic text diagnostics.

## Implementation structure

- Tests deterministic parsing, text, and safe reporting utilities. Syntax:
- result=test_dynamic_parse_text_reporting_suite()
- result -regression test result with explanatory messages
- The test checks operator-specification parsing, isotope predicates,
- label lookup, silent reporting calls, and polyadic text diagnostics.
- Announce the test target
- State the utility target of the test
- Build a small spin-system descriptor for parsing helpers
- Check isotope and spin-label parsing into single-spin opspecs
- Check product-operator parsing and Lx expansion coefficients
- Check label lookup against a unique label list
- Check electron and nucleus isotope predicates
