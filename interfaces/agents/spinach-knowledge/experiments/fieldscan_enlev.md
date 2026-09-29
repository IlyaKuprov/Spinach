# experiments/fieldscan_enlev.m

- MATLAB implementation: [experiments/fieldscan_enlev.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/fieldscan_enlev.m)

Source: https://spindynamics.org/wiki/index.php?title=fieldscan_enlev.m

Signature: fieldscan_enlev(spin_system,parameters). The routine makes a field-dependent energy-level plot; it does not return a numerical array or a magnetisation signal.

At the fixed requested orientation it obtains lab-frame Zeeman and coupling Hamiltonians and forms H(B)=B*Hz+Hc for each point of a linearly spaced field grid. At each field, Arnoldi/eigs is used to obtain the requested low-energy eigenvalues; the real parts are sorted and converted from angular-frequency units to inverse centimetres before plotting energy against field. This is a level diagram that can help inspect field-dependent level mixing or crossings, not a DNP or hyperpolarisation simulation and not an imaging calculation.

Required inputs:

- parameters.fields: two ascending endpoints [from,to] in tesla; the generated grid includes both endpoints.
- parameters.npoints: number of field-grid points (positive integer).
- parameters.orientation: fixed orientation as three Euler angles [alp bet gam] in radians.
- parameters.nstates: number of lowest energy levels to calculate (positive integer).
- spin_system: must use zeeman-hilb formalism and be built with `sys.magnet=1` T; the grumbler rejects any other stored reference field before the scan.

Output: a figure with magnetic field in tesla on the horizontal axis and energy in cm^-1 on the vertical axis. There is no returned signal vector. The source gives no fixed numerical field range, point count, orientation, or state-count example; those are user-selected inputs.

Limit: the plot is based on the specified spin Hamiltonian and fixed orientation. It does not propagate a density operator, calculate populations, include a sweep-rate response, or model DNP polarisation transfer.
