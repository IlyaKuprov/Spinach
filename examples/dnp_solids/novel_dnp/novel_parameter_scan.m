% 2D parameter scan of a NOVEL DNP experiment. <I_z> on 1H after a 
% 0.25 us contact time is calculated as a function of electron pul-
% se amplitude and offset. Further information in:
%
%                https://doi.org/10.1063/1.5000528
%
% Calculation time: minutes
%
% venkata-subbarao.redrouthu@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de

function novel_parameter_scan()

% X-band magnet
sys.magnet=0.34;

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
parameters.timestep=1e-9;                % Sequence time step, seconds
parameters.nsteps=250;                   % Number of time steps
parameters.flippulse=1;                  % NOVEL pulse sequence
parameters.needs={'aniso_eq'};           % Sequence needs rho_eq
parameters.grid='rep_2ang_100pts_sph';   % Spherical averaging grid

% Offset range, Hz
offsets=linspace(-35e6,35e6,71);
reference_point=-3.3e6;

% Nutation frequency range, Hz
nutfrqs=linspace(1e6,30e6,30);

% Hush the output
spin_system.sys.output='hush';

% Generate axis ticks
[X,Y]=meshgrid(offsets/1e6,nutfrqs/1e6);

% Get the figure going
figure(); dnp_surf=zeros(numel(nutfrqs),numel(offsets));

% Loop over powers
for n=1:numel(nutfrqs)

    % Set nutation frequency
    parameters.irr_powers=nutfrqs(n);

    % Loop over offsets
    parfor m=1:numel(offsets)

        % Localise parameters
        localpar=parameters;

        % Set the offsets
        localpar.offset=[offsets(m)+reference_point 0];

        % Run the simulation
        contact_curve=powder(spin_system,@noveldnp,localpar,'esr');
        dnp_surf(n,m)=real(contact_curve(end));

    end

    % Update the plot
    contourf(X,Y,dnp_surf,100,'edgecolor','none');
    kylabel('Electron nutation frequency, MHz');
    kxlabel('Microwave resonance offset, MHz');
    kcolourbar('$I_\textrm{z}$ expectation value on $^{1}$H');
    shading flat; colormap jet; drawnow;

end 

end

