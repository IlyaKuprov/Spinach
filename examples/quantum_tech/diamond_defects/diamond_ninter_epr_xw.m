% Field-swept powder EPR spectra of a nitrogen interstitial
% defect in diamond at X and W bands.
%
% alexey.bogdanov@weizmann.ac.il

function diamond_ninter_epr_xw()

% Set WAR9 nitrogen interstitial centre parameters
ninter_params.orientation='111';
ninter_params.nitrogen='15N';
ninter_params.centre='war9';

% Build the spin system
[sys,inter]=diamond_n_inter(ninter_params);

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
parameters.grid=6;
parameters.fwhm=1e-5;
parameters.int_tol=0.01;
parameters.tm_tol=0.01;
parameters.npoints=1024;
parameters.rspt_order=Inf;

% Set X-band parameters
parameters.mw_freq=9.755e9;
parameters.window=[0.347 0.349];

% Run the X-band simulation
[spec_x,par_x]=fieldsweep(spin_system,parameters);

% Plot the X-band spectrum
kfigure(); scale_figure([1.50 0.75]);
subplot(1,2,1); plot(par_x.b_axis,spec_x);
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
ktitle('WAR9 N interstitial X-band EPR');
xlim tight; ylim padded; kgrid;

% Set W-band parameters
parameters.mw_freq=94e9;
parameters.window=[3.351 3.355];

% Run the W-band simulation
[spec_w,par_w]=fieldsweep(spin_system,parameters);

% Plot the W-band spectrum
subplot(1,2,2); plot(par_w.b_axis,spec_w);
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
ktitle('WAR9 N interstitial W-band EPR');
xlim tight; ylim padded; kgrid;

end

