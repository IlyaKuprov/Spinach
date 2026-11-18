# gcce_finite_drive

Syntax: [answer,diagnostics]=gcce_finite_drive().

Small electron–carbon model demonstrating a 20-microsecond rectangular observed
spin drive at 40000 rad/s followed by free evolution. Canonical sys/inter,
transition preparation, a scalar parameters.drive descriptor, and generalised
clexpn with @finite_drive avoid manually assembled states and event arrays.
The returned answer is the column-vectorised observed density, one row per time;
diagnostics contains reduction and error information. The plot shows populations.
The independent two-spin expm reference lives in test_cce_examples.
