# tests/kernel/test_dynamic_overload_ttclass_suite.m

- Signature: `result=test_dynamic_overload_ttclass_suite()`

Regression test comparing `ttclass` overloads with dense references for deterministic one- and two-core tensor trains. Covers construction, shape and indexing, arithmetic, multiplication and inner products, complex operations, reductions, and vectorisation. Also checks packing, orthogonalisation, truncation, shrinkage, AMEn summation and scalar solving, `save_anyway` round-trip, random-rank construction, and rejection of a non-tensor right-hand side by `mldivide`. Returns a test result with explanatory messages.
