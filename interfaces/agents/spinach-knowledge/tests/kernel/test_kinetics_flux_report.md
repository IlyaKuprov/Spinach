# tests/kernel/test_kinetics_flux_report.m

Tests the original directed point-to-point flux table and flux-type diagnostics in both `create` and explicit `kinetics(...,'report')` output. Asymmetric rates expose reversed from/to indices; intramolecular and intermolecular cases each require exactly one table per call. Report mode must leave the chemical kinetics generator unchanged.

Absent and all-zero flux inputs must remain quiet, with identical ordinary/report-mode chemical generators. This is a reporting regression, not a replacement for the flux-generator tests.
