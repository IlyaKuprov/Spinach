% The transformation of -E_z into I_z during the contact time of the
% X-inverse-X DNP experiment. Further information in:
%
%               https://doi.org/10.1021/jacs.1c09900
%
% Calculation time: seconds
%
% venkata-subbarao.redrouthu@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de

function xix_contact_curve()

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
parameters.offset=[(-13.0+59.5)*1e6 0];  % -13 MHz reference point, 59.5 MHz offset
parameters.irr_powers=17.8e6;            % Electron nutation frequency [Hz]
parameters.pulse_dur=48e-9;              % Pulse duration, seconds
parameters.nloops=80;                    % Number of XiX DNP blocks
parameters.phase=pi;                     % Second pulse inverted phase
parameters.grid='rep_2ang_1600pts_sph';  % Powder averaging grid
parameters.needs={'aniso_eq'};           % Sequence needs rho_eq

% Run the calculation
contact_curve=powder(spin_system,@xixdnp,parameters,'esr');

% Plotting
time_axis=linspace(0,2*parameters.pulse_dur*parameters.nloops,parameters.nloops+1);
figure(); plot(time_axis,real(contact_curve)); kxlabel('contact time, seconds');
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H'); xlim tight; kgrid;

end

