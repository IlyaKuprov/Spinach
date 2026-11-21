# tests/kernel/test_cce_no_legacy.m

`result=test_cce_no_legacy()` checks that the retired CCE driver, wrappers, and uncalled raw-control interpreter cannot resolve as functions on the active Spinach path. MATLAB directory status 7 is allowed because kernel/cce remains a directory, but which must still return empty. It locates the checkout from its own filename and scans Git-tracked MATLAB files, excluding itself, for direct calls and function handles after stripping comments and quoted text. Transpose operators remain executable text.

The test also requires the native context, bath constructor, orientation composition, and registered experiments to resolve inside the same checkout. Scanner controls distinguish inert quoted references from a call following a transpose. No physical signal or numerical tolerance is changed.

The test supplies all explicitly removed inventory paths to clexpn, together with plausible relocated spellings, and requires native option errors. It also tests all eight direct experiments against finite, ideal, empty-structure, and empty-numeric raw event arrays, requiring an error naming native controls.
