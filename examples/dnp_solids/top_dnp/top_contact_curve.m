% The transformation of -E_z into I_z during the contact time of the
% time-optimised pulsed DNP experiment. Further information in:
%
%              https://doi.org/10.1126/sciadv.aav6909
%
% Calculation time: seconds
%
% venkata-subbarao.redrouthu@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de

function top_contact_curve()

% Q-band magnet
sys.magnet=1.2142;

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
parameters.offset=[(-13.0+92.5)*1e6 0];  % -13 MHz reference point, 92.5 MHz offset
parameters.irr_powers=17.8e6;            % electron nutation frequency [Hz]
parameters.pulse_dur=10e-9;              % pulse duration, seconds
parameters.delay_dur=14e-9;              % delay duration, seconds
parameters.nloops=300;                   % number of TOP DNP blocks
parameters.grid='rep_2ang_3200pts_sph';  % Spherical averaging grid
parameters.needs={'aniso_eq'};           % Sequence needs rho_eq

% Run the calculation
contact_curve=powder(spin_system,@topdnp,parameters,'esr');

% Plotting
time_axis=linspace(0,parameters.nloops*...
                    (parameters.pulse_dur+parameters.delay_dur),...
                     parameters.nloops+1);
figure(); plot(time_axis,real(contact_curve)); kxlabel('contact time, seconds');
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H'); axis tight; kgrid;

end

