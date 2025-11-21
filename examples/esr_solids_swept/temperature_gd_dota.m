% Powder averaged W-band field-swept ESR spectrum of Gd(III)
% DOTA complex. Exact diagonalisation is used and a tempera-
% ture dependence plot is produced.
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il

function temperature_gd_dota()

% Isotopes
sys.isotopes={'E8'};

% Magnet field (must be 1)
sys.magnet=1;

% Properties
inter.zeeman.scalar={1.9918};
inter.coupling.eigs{1,1}=[0.57e9 0.57e9 -2*0.57e9]/3;
inter.coupling.euler{1,1}=[0 0 0];

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Get figure going
kfigure(); hold on; n=1;
scale_figure([2.0 1.5]);

% Loop over temperatures
for T=[100 10 1 0.1]

    % Set the temperature
    inter.temperature=T;

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Experiment parameters
    parameters.spins={'E8'};
    parameters.grid=4;
    parameters.mw_freq=90e9;
    parameters.fwhm=2e-4;
    parameters.int_tol=0.001;
    parameters.tm_tol=0.1;
    parameters.window=[3.05 3.4];
    parameters.npoints=4096;
    parameters.rspt_order=Inf;

    % Run the simulation
    [b_axis,spec]=fieldsweep(spin_system,parameters);

    % Plotting
    subplot(2,2,n); plot(b_axis',spec');
    kxlabel('magnetic field, tesla');
    kylabel('intensity, a.u.'); kgrid;
    ktitle([num2str(T) ' K']);
    xlim tight; drawnow; n=n+1;

end

end

