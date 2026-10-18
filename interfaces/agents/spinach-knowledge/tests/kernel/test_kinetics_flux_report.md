# tests/kernel/test_kinetics_flux_report.m

Tests creation-time reporting of asymmetric directed point-to-point flux rates for intramolecular and intermolecular transfer, including source/destination orientation. Absent and zero flux produce no flux table. After basis construction, `kinetics(spin_system)` is quiet and repeated calls return identical generators.
