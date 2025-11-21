% PSYCOSY of one salsalate ring.
%
% Calculation time: minutes on NVidia Tesla A100, much longer on CPU
% 
% a.m.kenwright@durham.ac.uk

function psycosy_salsalate()

% Magnet
sys.magnet=14.1;

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H'};

% Interactions
inter.zeeman.scalar={8.14 7.44 7.71 7.28 7.60};
inter.coupling.scalar{1,2}=7.9; 
inter.coupling.scalar{2,3}=7.5; 
inter.coupling.scalar{3,4}=8.1; 
inter.coupling.scalar{1,3}=1.6; 
inter.coupling.scalar{4,2}=1.2;
inter.coupling.scalar{5,2}=0;
inter.coupling.scalar{5,5}=0;

% Algorithmic options
sys.enable={'greedy'}; % 'gpu'
sys.tols.prox_cutoff=4.0;
sys.tols.merge_dim=500;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.space_level=1;
bas.connectivity='scalar_couplings';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sample geometry
parameters.dims=0.015;    % m
parameters.npts=100;
parameters.deriv={'period',3};

% Diffusion and flow
parameters.u=zeros(parameters.npts,1);
parameters.diff=0;

% Relaxation phantom
parameters.rlx_ph={zeros(parameters.npts,1)};
parameters.rlx_op={relaxation(spin_system)};

% Initial and detection state phantoms
parameters.rho0_ph={ones(parameters.npts,1)};
parameters.rho0_st={state(spin_system,'Lz','1H')};
parameters.coil_ph={ones(parameters.npts,1)};
parameters.coil_st={state(spin_system,'L+','1H')};

% Sequence parameters
parameters.offset=4620;
parameters.sweep=720;
parameters.npoints=[512 512];
parameters.zerofill=[1024 1024];
parameters.spins={'1H'};
parameters.axis_units='ppm';
parameters.tmix=110e-3;       % mixing time, seconds
parameters.gamp=1e-2;         % gradient amplitude, T/m

% Saltire chirp parameters
parameters.sal_ang=20;        % flip angle of the saltire chirp (degrees)
parameters.sal_dur=0.015;     % pulse width of saltire chirp (s)
parameters.sal_del=0.05;      % chirp pulse gradient duration (s)
parameters.sal_swp=10000;     % sweep width of saltire chirp (Hz)
parameters.sal_npt=250;       % number of points in the saltire chirp
parameters.sal_smf=20;        % saltire chirp smoothing factor

% Simulation
fid=imaging(spin_system,@psycosy,parameters);

% Apodisation
fid=apodisation(spin_system,fid,{{'sqsin'},{'sqsin'}});

% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill(2),parameters.zerofill(1)));

% Plotting 
kfigure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,abs(spectrum),parameters,...
        20,[0.02 0.2 0.02 0.2],2,256,6,'positive');
              
end

