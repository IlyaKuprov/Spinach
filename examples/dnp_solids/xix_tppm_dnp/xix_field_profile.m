% Field profile of a XiX DNP experiment. <I_z> after a fixed 
% contact time is calculated as a function of electron pulse
% amplitude and offset. Further information in:
%
%            https://doi.org/10.1021/jacs.1c09900
%
% Calculation time: minutes (a large powder grid is needed)
%
% venkata-subbarao.redrouthu@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de

function xix_field_profile()

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

% Hush the output
sys.output='hush';

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
parameters.irr_powers=17.8e6;           % electron nutation frequency [Hz]
parameters.pulse_dur=48e-9;             % pulse duration, seconds
parameters.nloops=150;                  % number of XiX DNP blocks
parameters.phase=pi;                    % Second pulse in opposite phase
parameters.grid='rep_2ang_1600pts_sph'; % Spherical averaging grid
parameters.needs={'aniso_eq'};          % Sequence needs rho_eq

% MW resonance offset grid, Hz
offsets=linspace(-150e6,150e6,120);
reference_point=-13e6;

% Parallel loop over offsets
parfor m=1:numel(offsets)

    % Localise parameters array
    localpar=parameters;

    % Set the offsets
    localpar.offset=[offsets(m)+reference_point 0];

    % Get the contact curve
    contact_curve=powder(spin_system,@xixdnp,localpar,'esr');
    
    % Record the last point
    field_prof(m)=real(contact_curve(end));

end 

% Plotting
figure(); plot(offsets/1e6,field_prof); axis tight;
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');  
kxlabel('Microwave resonance offset, MHz'); kgrid;

end

