% Field-swept powder EPR spectra of Ni defects
% in diamond at X and W bands.
%
% alexey.bogdanov@weizmann.ac.il

function diamond_ni_epr_xw()

% Set Ni NE1 centre model parameters.
ni_params.centre='ne1';
ni_params.orientation='111';
ni_params.nickel='61Ni';

% Build the spin system
[sys,inter]=diamond_ni(ni_params);

% Field sweep
sys.magnet=1;

% Define the basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% EPR sim parameters
parameters.spins={'E'};
parameters.grid=20;
parameters.fwhm=0.0001;
parameters.int_tol=0.1;
parameters.tm_tol=0.01;
parameters.npoints=2048;
parameters.rspt_order=Inf;

% Set X-band parameters
parameters.mw_freq=9.5e9;
parameters.window=[0.3 0.36];

% Run the X-band simulation
[b_axis_x,spec_x]=fieldsweep(spin_system,parameters);

% Plot the X-band spectrum
kfigure(); scale_figure([1.50 0.75]);
subplot(1,2,1); plot(b_axis_x',spec_x');
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
ktitle('Ni NE1 X-band EPR');
xlim tight; ylim padded; kgrid;

% Set W-band parameters
parameters.mw_freq=94e9;
parameters.window=[3.1 3.4];

% Run the W-band simulation
[b_axis_w,spec_w]=fieldsweep(spin_system,parameters);

% Plot the W-band spectrum
subplot(1,2,2); plot(b_axis_w',spec_w');
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
ktitle('Ni NE1 W-band EPR');
xlim tight; ylim padded; kgrid;

end

