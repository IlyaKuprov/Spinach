% 4Q ultrafast MaxQ NMR spectrum for a coupled four-spin 
% system in the presence of realistic diffusion.
%
% Calculation time: minutes on NVidia Tesla A100, much longer on CPU
%
% mariagrazia.concilio@sjtu.edu.cn
% ilya.kuprov@weizmann.ac.il
% jean-nicolas.dumez@univ-nantes.fr

function ufmq_4spin()

% Magnetic field
sys.magnet=14.1; % Tesla

% Chemical shifts
sys.isotopes={'1H','1H','1H','1H'};
inter.zeeman.scalar={0.50,0.35,0.15,0};
 
% 3J couplings  
inter.coupling.scalar{1,2}=8.0;
inter.coupling.scalar{2,3}=8.0;
inter.coupling.scalar{3,4}=8.0;
  
% 4J couplings  
inter.coupling.scalar{1,3}=3.0;
inter.coupling.scalar{2,4}=3.0;
 
% 5J couplings  
inter.coupling.scalar{1,4}=2.0;
inter.coupling.scalar{4,4}=0;
 
% Coherence selection
parameters.mqorder=+4;

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
parameters.dims=0.015;         % m
parameters.npts=500;
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
parameters.diff=18e-10;      % m^2/s

% Acquisition parameters
parameters.spins={'1H'};
parameters.offset=0;
parameters.deltat=6e-06;  % s
parameters.npoints=120;
parameters.nloops=50;

% Calculate maximal k-value
k_max=parameters.npoints/parameters.dims;

% Get acquisition gradient duration
Ta=parameters.deltat*parameters.npoints; % s

% Get acquisition gradient amplitude
parameters.Ga=((2*pi)*k_max)/(spin(parameters.spins{1})*Ta); % T/m

% Timings
parameters.delay=0.041;   % s

% Encoding parameters
parameters.pulsenpoints=500; 
parameters.nWURST=40;
parameters.Te=0.015;   % s
parameters.BW=15000;   % Hz
parameters.Ge=0.023;   % T/m
parameters.chirptype='wurst';

% Simulation
ktdata=imaging(spin_system,@ufmq,parameters);

% Plot echoes
kfigure(); subplot(1,2,1);
imagesc(imag(ktdata));
xlabel('t_{2} dimension points'); 
ylabel('k-space points');

% Fourier transform in the conventional dimension
kwdata_sim=fftshift(fft(ktdata,[],2),2);

% Plot the spectrum
subplot(1,2,2);
parameters.axis_units='ppm';
parameters.offset_uf_cov=0;
parameters.spins={'1H','1H'};
parameters.offset=[0 0];
plot_uf(spin_system,abs(kwdata_sim),parameters);

end

