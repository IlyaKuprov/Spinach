% Field-swept powder EPR spectra of a P1 centre 
% in diamond at X and W bands.
%
% alexey.bogdanov@weizmann.ac.il

function diamond_p1_epr_xw()

% Set P1 model parameters.
p1_params.orientation='111';
p1_params.nitrogen='14N';

% Build the spin system.
[sys,inter]=diamond_p1(p1_params);

% Field sweep
sys.magnet=1;

% Define the basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set common EPR parameters
parameters.spins={'E'};
parameters.grid=6;
parameters.fwhm=1e-4;
parameters.int_tol=0.01;
parameters.tm_tol=0.1;
parameters.npoints=1024;
parameters.rspt_order=Inf;

% Set X-band parameters
parameters.mw_freq=9.5e9;
parameters.window=[0.33 0.35];

% Run the X-band simulation
[b_axis_x,spec_x]=fieldsweep(spin_system,parameters);

% Plot the X-band spectrum
kfigure(); scale_figure([1.50 0.75]);
subplot(1,2,1); plot(b_axis_x',spec_x');
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
ktitle('P1 X-band EPR');
xlim tight; ylim padded; kgrid;

% Set W-band parameters
parameters.mw_freq=94e9;
parameters.window=[3.348 3.36];

% Run the W-band simulation
[b_axis_w,spec_w]=fieldsweep(spin_system,parameters);

% Plot the W-band spectrum
subplot(1,2,2); plot(b_axis_w',spec_w');
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
ktitle('P1 W-band EPR');
xlim tight; ylim padded; kgrid;

end
