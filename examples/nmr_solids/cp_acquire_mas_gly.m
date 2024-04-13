% 1H-13C cross-polarisation followed by acquisition under magic
% angle spinning in alpha-glycine powder. Reduced Liouville spa-
% ce is used: up to, and including, three-spin correlations.
%
% Calculation time: minutes on Tesla A100, much longer on CPU.
%
% guinevere.mathies@uni-konstanz.de

function cp_acquire_mas_gly()

% Spin system properties (PCM DFT calculation)
[sys,inter]=g2spinach(gparse('..\..\examples\standard_systems\glycine.log'),...
                      {{'H','1H'},{'C','13C'}},[31.8 182.1],[]);

% 400 MHz spectrometer
sys.magnet=9.4;

% Isotropic alpha-glycine chemical shifts
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1,176.4); % CO
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,2, 43.6); % CA
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,3,  2.6); % H_CA
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,4,  3.8); % H_CA
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,5,  8.0); % H_N
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,6,  8.0); % H_N
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,7,  8.0); % H_N

% Spin temperature
inter.temperature=298;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=3;

% Algorithmic options
sys.enable={'greedy','gpu'};
sys.disable={'pt'};

% Neglect interactions below 200 Hz
sys.tols.inter_cutoff=2*pi*200;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.spins={'1H','13C'};
parameters.rate=10000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=5;
parameters.grid='rep_2ang_100pts_sph';
parameters.offset=[2e3 10e3];  % Hz
parameters.hi_pwr=83e3;        % Hz
parameters.cp_pwr=[60e3 50e3]; % Hz
parameters.cp_dur=50e-5;       % seconds
parameters.sweep=5e4;          % Hz
parameters.npoints=512;
parameters.zerofill=4096;
parameters.needs={'iso_eq'};   % initial condition
parameters.verbose=1;

% Detection state
parameters.coil=state(spin_system,'L+','13C');

% Simulation
fid=singlerot(spin_system,@cp_acquire_soft,parameters,'nmr');

% Apodization
fid=apodization(fid,'exp-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% 1D plotting
parameters.offset=parameters.offset(2);
parameters.spins=parameters.spins(2); figure();
plot_1d(spin_system,real(spectrum),parameters);

end

