% X-band pulse-acquire FFT ESR spectrum of a bisnitroxide radical, 
% using explicit time domain simulation with Redfield relaxation
% supeoperator. Parameters from https://doi.org/10.1039/C8CP06819D 
%
% Calculation time: seconds
%
% fmentink@magnet.fsu.edu

function relaxation_bisnitroxide()

% Spin system
sys.isotopes={'E','E','14N','14N'};

% Zeeman interactions
inter.zeeman.eigs{1}=[2.00925 2.00605 2.00205];
inter.zeeman.euler{1}=[0 0 0];
inter.zeeman.eigs{2}=inter.zeeman.eigs{1};
inter.zeeman.euler{2}=[123.1 129.8 -46.6]*pi/180;
inter.zeeman.eigs{3}=[0 0 0];
inter.zeeman.euler{3}=[0 0 0];
inter.zeeman.eigs{4}=[0 0 0];
inter.zeeman.euler{4}=[0 0 0];

% Couplings
inter.coupling.eigs=cell(4,4);
inter.coupling.euler=cell(4,4);
inter.coupling.scalar=cell(4,4);
inter.coupling.eigs{1,2}=[17.5,17.5,-35]*1e6;
inter.coupling.eigs{1,3}=[18,17,103]*1e6;
inter.coupling.eigs{2,4}=inter.coupling.eigs{1,3};
inter.coupling.euler{1,2}=[-174 74 0]*pi/180;
inter.coupling.euler{1,3}=[0 0 0]*pi/180;
inter.coupling.euler{2,4}=inter.zeeman.euler{2};
inter.coupling.scalar{1,2}=-2*16e6;

% Magnet induction
sys.magnet=0.35; % Telsa

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={4e-10};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E');
parameters.coil=state(spin_system,'L+','E');
parameters.decouple={};
parameters.offset=0e8;
parameters.sweep=5e8;
parameters.npoints=512;
parameters.zerofill=1024;
parameters.axis_units='GHz-labframe';
parameters.derivative=1;
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'esr');

% Apodization
fid=apodization(fid,'none-1d');

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

