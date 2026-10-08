% The transformation of -E_z into I_z during the contact time of the
% BEAM DNP experiment. Further information in:
%
%               https://doi.org/10.1126/sciadv.abq0536
%
% Calculation time: seconds
%
% venkata-subbarao.redrouthu@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de

function beam_contact_curve()

% X-band magnet
sys.magnet=0.3483;

% Electron and two protons
sys.isotopes={'E','1H','1H'};

% Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
inter.zeeman.eigs={[2.00319 2.00319 2.00258],[0 0 5],[0 5 0]};
inter.zeeman.euler=(pi/180)*{[0 10 0],[0 0 10],[100 0 0]};

% Cartesian coordinates
inter.coordinates={[0.000 0.000 0.000];
                   [0.000 3.500 0.000];
                   [2.475 2.475 0.000]};

% Spin temperature
inter.temperature=80;

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Detection state
parameters.coil=state(spin_system,'Lz','1H');

% Experiment parameters
parameters.spins={'E','1H'};
parameters.offset=[(-3.3+5.0)*1e6 0];    % -13 MHz reference point, 5.0 MHz offset
parameters.irr_powers=32.0e6;            % Electron nutation frequency [Hz]
parameters.pulse_dur=[20.0e-9 28.7e-9];  % Pulse durations, seconds
parameters.nloops=165;                   % Number of BEAM DNP blocks
parameters.grid='rep_2ang_800pts_sph';   % Powder averaging grid
parameters.needs={'aniso_eq'};           % Sequence needs rho_eq

% Run the calculation
contact_curve=powder(spin_system,@beamdnp,parameters,'esr');

% Plotting
time_axis=linspace(0,sum(parameters.pulse_dur)*parameters.nloops,parameters.nloops+1);
kfigure(); plot(time_axis,real(contact_curve)); kxlabel('contact time, seconds');
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H'); xlim tight; kgrid;

end

