# tests/kernel/test_kinetics_flux_report.m

Tests creation-time reporting of asymmetric directed point-to-point flux rates for intramolecular and intermolecular transfer, including source/destination orientation. Absent and zero flux produce no flux table. After basis construction, `kinetics(spin_system)` is quiet and repeated calls return identical generators.

Current interface note: flux rates/type live at inter.flux_rate and inter.flux_type, then spin_system.rlx.flux_rate and spin_system.rlx.flux_type. Selecting 'flux' in inter.relaxation activates the contribution; tests of prior automatic or chemistry-owned behavior are historical, not current API assertions.
