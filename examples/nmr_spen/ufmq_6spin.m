% 6Q ultrafast MaxQ NMR spectrum for a coupled six-spin 
% system in the presence of realistic diffusion.
%
% Calculation time: hours on NVidia Tesla A100, much longer on CPU
%
% mariagrazia.concilio@sjtu.edu.cn
% ilya.kuprov@weizmann.ac.il
% jean-nicolas.dumez@univ-nantes.fr

function ufmq_6spin()

% Magnetic field
sys.magnet=14.1; % Tesla

% Chemical shifts
sys.isotopes={'1H','1H','1H','1H','1H','1H'};
inter.zeeman.scalar={-1.0 -0.5 0.0 +0.3 +0.7 +0.9};

% 3J couplings  
inter.coupling.scalar{1,2}=8.00;
inter.coupling.scalar{2,3}=8.00;
inter.coupling.scalar{3,4}=8.00;
inter.coupling.scalar{4,5}=8.00;

% 4J couplings  
inter.coupling.scalar{1,3}=4.00;
inter.coupling.scalar{1,4}=4.00;
inter.coupling.scalar{2,4}=4.00;
inter.coupling.scalar{3,5}=4.00;
inter.coupling.scalar{1,6}=4.00;
inter.coupling.scalar{5,6}=4.00;

% 5J couplings 
inter.coupling.scalar{2,5}=4.00;
inter.coupling.scalar{2,6}=4.00;
inter.coupling.scalar{4,6}=4.00;

% 6J couplings
inter.coupling.scalar{1,5}=2.00;
inter.coupling.scalar{3,6}=2.00;
inter.coupling.scalar{6,6}=0;
 
% Coherence selection
parameters.mqorder=+6;

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
parameters.npts=300;
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
figure(2); subplot(1,2,1);
imagesc(real(ktdata));
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

