% Field profile of a BEAM DNP experiment. <I_z> after a fixed 
% contact time is calculated as a function of electron pulse
% amplitude and offset. Further information in:
%
%            https://doi.org/10.1126/sciadv.abq0536
%
% Calculation time: minutes (a large powder grid is needed)
%
% venkata-subbarao.redrouthu@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de

function beam_field_profile()

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
parameters.irr_powers=32.0e6;            % Electron nutation frequency [Hz]
parameters.pulse_dur=[20.0e-9 28.7e-9];  % Pulse durations, seconds
parameters.nloops=165;                   % Number of BEAM DNP blocks
parameters.grid='rep_2ang_800pts_sph';   % Powder averaging grid
parameters.needs={'aniso_eq'};           % Sequence needs rho_eq

% MW resonance offset grid, Hz
offsets=linspace(-60e6,60e6,120);
reference_point=-3.3e6;

% Parallel loop over offsets
parfor m=1:numel(offsets)

    % Localise parameters array
    localpar=parameters;

    % Set the offsets
    localpar.offset=[offsets(m)+reference_point 0];

    % Get the contact curve
    contact_curve=powder(spin_system,@beamdnp,localpar,'esr');
    
    % Record the last point
    field_prof(m)=real(contact_curve(end));

end 

% Plotting
kfigure(); plot(offsets/1e6,field_prof); axis tight;
kylabel('$I_\textrm{z}$ expectation value on $^{1}$H');  
kxlabel('Microwave resonance offset, MHz'); kgrid;

end

