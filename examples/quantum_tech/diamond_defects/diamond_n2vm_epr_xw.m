% Field-swept powder EPR spectra of an N2V- centre
% in diamond at X and W bands.
%
% alexey.bogdanov@weizmann.ac.il

function diamond_n2vm_epr_xw()

% Set N2V- centre parameters
n2v_params.orientation='111';
n2v_params.nitrogen='15N';
n2v_params.include_13c=false;

% Build the spin system
[sys,inter]=diamond_n2vm(n2v_params);

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
parameters.fwhm=0.00003;
parameters.int_tol=0.001;
parameters.tm_tol=0.01;
parameters.npoints=1024;
parameters.rspt_order=Inf;

% Set X-band parameters
parameters.mw_freq=9.755e9;
parameters.window=[0.347 0.349];

% Run the X-band simulation
[b_axis_x,spec_x]=fieldsweep(spin_system,parameters);


% Plot the X-band spectrum
kfigure(); scale_figure([1.50 0.75]);
subplot(1,2,1); plot(b_axis_x,spec_x);
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
ktitle('N2V$^{-}$ X-band EPR');
xlim tight; ylim padded; kgrid;

% Set W-band parameters
parameters.mw_freq=94e9;
parameters.window=[3.351 3.355];

% Run the W-band simulation
[b_axis_w,spec_w]=fieldsweep(spin_system,parameters);

% Plot the W-band spectrum
subplot(1,2,2); plot(b_axis_w,spec_w);
kxlabel('magnetic field, tesla');
kylabel('intensity, a.u.');
ktitle('N2V$^{-}$ W-band EPR');
xlim tight; ylim padded; kgrid;

end

