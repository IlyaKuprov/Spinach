% PSYCHE pure-shift NMR spectrum of strychnine.
%
% Calculation time: hours on NVidia Tesla A100, much longer on CPU
%
% mohammadali.foroozandeh@chem.ox.ac.uk
% ilya.kuprov@weizmann.ac.il

function psyche_strychnine()

% Strychnine spin system
[sys,inter]=strychnine({'1H'});

% Move shifts into spectral window
inter.zeeman.scalar=inter.zeeman.scalar-5.0;

% Magnetic induction
sys.magnet=11.7;

% Basis set
bas.formalism='sphten-liouv';
bas.connectivity='scalar_couplings';
bas.approximation='IK-2';
bas.space_level=1;

% Algorithmic options
sys.tols.inter_cutoff=2.0;
sys.disable={'pt','colorbar'};
sys.enable={'greedy'}; % 'gpu'

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.offset=0;
parameters.sweep=[100 5000]; % Hz
parameters.npoints=[32 2048];
parameters.zerofill=[128 8192];
parameters.spins={'1H'};
parameters.axis_units='Hz';
parameters.g_amp=0.015;  % T/m

% Saltire chirp parameters
parameters.beta=20;             % flip angle of the saltire chirp (degrees)
parameters.duration=0.015;      % pulse width of saltire chirp (s)
parameters.del=0;               % diffusion delay (s)
parameters.bandwidth=10000;     % sweep width of saltire chirp (Hz)
parameters.pulsenpoints=1000;   % number of points in the saltire chirp
parameters.smfactor=20;         % smoothing factor
parameters.chirptype='smoothed'; 

% Coherent evolution timesteps
parameters.timestep1=1/parameters.sweep(1);
parameters.timestep2=1/parameters.sweep(2);
parameters.delta=parameters.timestep1/4;

% Sample parameters
parameters.dims=0.015;  % m
parameters.npts=1000;
parameters.deriv={'period',7};

% Diffusion and flow
parameters.diff=0;      % m^2/s
parameters.u=zeros(parameters.npts,1);

% Initial and detection state phantoms
parameters.rho0_ph={ones(parameters.npts,1)};
parameters.rho0_st={state(spin_system,'Lz','1H','cheap')};
parameters.coil_ph={ones(parameters.npts,1)};
parameters.coil_st={state(spin_system,'L+','1H','cheap')};

% Relaxation phantom
parameters.rlx_ph={};
parameters.rlx_op={};

% Simulation
fid=imaging(spin_system,@psyche,parameters);

% Reconstructing the pure shift FID
np_chunk=parameters.sweep(2)/parameters.sweep(1);
fidps=fid(1:np_chunk,:); fidps=fidps(:);

% Apodisation
fidps=apodisation(spin_system,fidps,{{'gauss',6}});

% Fourier transform
spectrum_2d=fftshift(fft2(fid,parameters.zerofill(2),...
                              parameters.zerofill(1)));
spectrum_ps=fftshift(fft(fidps,parameters.zerofill(2)));

% Plotting: full 2D version
figure(); scale_figure([2.0 1.5]); subplot(1,2,1);
plot_2d(spin_system,abs(spectrum_2d),parameters,...
        20,[0.05 1.0 0.05 1.0],2,256,6,'positive');
 
% Plotting: 1D projection
subplot(1,2,2);
parameters.sweep=parameters.sweep(2);
parameters.zerofill=parameters.zerofill(2);
plot_1d(spin_system,imag(spectrum_ps),parameters);

end
