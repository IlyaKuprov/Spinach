# tests/kernel/test_dynamic_parse_text_reporting_suite.m

- Signature: `result=test_dynamic_parse_text_reporting_suite()`

Regression checks for operator-specification parsing, isotope predicates, label lookup, and text reporting. The suite exercises `human2opspec` selection and product-operator coefficients, `idxof` label lookup, and electron/nucleus classification. It checks that `report`, `banner`, and `summary_coordinates` are silent in hush mode, and that `polinfo` reports the shape of a polyadic product and its matrix-core sizes. These are described as test checks, not as verified passing results.
