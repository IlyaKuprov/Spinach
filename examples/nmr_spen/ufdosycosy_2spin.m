% Ultrafast 3D DOSY-COSY for a two-spin system.
%
% Calculation time: hours on NVidia Tesla A100, much longer on CPU
%
% Ludmilla Guduff
% Jean-Nicolas Dumez

function ufdosycosy_2spin()

% Spin system
sys.isotopes={'1H','1H'};

% Interactions
sys.magnet=14.1;
inter.zeeman.scalar={2.0,-1.3};
inter.coupling.scalar{1,2}=15;
inter.coupling.scalar{2,2}=0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.disable={'pt'};
sys.enable={'greedy','gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Assumptions
spin_system=assume(spin_system,'nmr');

% Sample geometry
parameters.dims=0.015;    % m
parameters.npts=3000;
parameters.deriv={'period',7};

% Relaxation phantom
parameters.rlx_ph={};
parameters.rlx_op={};

% Initial and detection state phantoms
parameters.rho0_ph={ones(parameters.npts,1)};
parameters.rho0_st={state(spin_system,'Lz','1H','cheap')};
parameters.coil_ph={ones(parameters.npts,1)};
parameters.coil_st={state(spin_system,'L+','1H','cheap')};

% Diffusion and flow
parameters.u=zeros(parameters.npts,1);
parameters.diff=4e-10;    % m^2/s

% Acquisition parameters
parameters.spins={'1H'};
parameters.td=0.100;      % s
parameters.deltat=1.5e-6; % s
parameters.sweep=1302;
parameters.npoints1=128;
parameters.npoints2=64;
parameters.nloops=256;
parameters.Ga=0.52;       % T/m 
parameters.offset=600;    % Hz

% Third dimension parameters 
parameters.sweep=1302;    % Hz

% Encoding parameters
parameters.pulsenpoints=1000; 
parameters.smfactor=0.1;
parameters.Te=0.0015;     % s
parameters.Tau=0.0016;    % s
parameters.BW=110000;     % Hz
parameters.Ge=0.2535;     % T/m 
parameters.chirptype='smoothed';

% Equation for the fit
parameters.dscale=1e-10;  % diffusion coefficient scale 
parameters.fovmin=-0.0041;
parameters.fovmax=+0.0041;

% Window multiplication and apodization 
parameters.apodize='hamming';
parameters.window='sine';

% Fitting model
parameters.model='keeler_corr';

% Simulation
fid=imaging(spin_system,@spendosycosy,parameters);

% Processing 
spectrum=fftshift(fft(fid,[],1),1);
spectrum=fftshift(fft(spectrum,[],2),2);
spectrum=fftshift(fft(spectrum,[],3),3);

% Plotting
figure(); volplot(abs(spectrum).^(1/2),[-1 1 -1 1 -1 1]);
    
end 

