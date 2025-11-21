% Ultrafast COSY for a coupled two-spin system.
%
% Calculation time: minutes on NVidia Tesla A100, much longer on CPU
%
% Jean-Nicolas Dumez
% Ludmilla Guduff

function ufcosy_2spin()

% Interactions
sys.magnet=14.0;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={3.70 4.50}-{4.0 4.0};
inter.coupling.scalar{1,2}=10;
inter.coupling.scalar{2,2}=0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.disable={'pt'};
sys.enable={'greedy'}; % 'gpu'

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sample geometry
parameters.dims=0.015;    % m
parameters.npts=500;
parameters.deriv={'period',3};

% Relaxation phantom
parameters.rlx_ph={zeros(parameters.npts,1)};
parameters.rlx_op={relaxation(spin_system)};

% Initial and detection state phantoms
parameters.rho0_ph={ones(parameters.npts,1)};
parameters.rho0_st={state(spin_system,'Lz','1H','cheap')};
parameters.coil_ph={ones(parameters.npts,1)};
parameters.coil_st={state(spin_system,'L+','1H','cheap')};

% Diffusion and flow
parameters.u=zeros(parameters.npts,1);
parameters.diff=0;        % m^2/s

% Acquisition parameters
parameters.spins={'1H'};
parameters.offset=0.0;
parameters.deltat=0.5e-6; % s
parameters.npoints=512;
parameters.nloops=128;
parameters.Ga=0.50;       % T/m

% Encoding parameters
parameters.pulsenpoints=1000; 
parameters.nWURST=40;
parameters.Te=0.015; % s
parameters.BW=10e3;  % Hz
parameters.Ge=0.01;  % T/m

% Coherence selection parameters
parameters.Gp=0.47;  % T/m
parameters.Tp=0.001; % s

% Simulation
fid=imaging(spin_system,@spencosy,parameters);

% Processing and plotting
spectrum=fftshift(fft(fid,[],2),2);
kfigure(); contour(abs(spectrum)); kgrid;

end

