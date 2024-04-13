% Two-site position exchange for a deuterium nucleus. The sites
% differ in the chemical shift and the orientation of the quad-
% rupolar tensor.
%
% Calculation time: seconds.
%
% umitakbey@gmail.com
% i.kuprov@soton.ac.uk

function mas_exchange_1()

% System specification
sys.magnet=9.4;

% Spin system
sys.isotopes={'2H','2H'};

% Quadrupolar interactions
[Q1,Q2]=weblab2nqi(0.16e6,0.1,1,0,acos(1/3),2*pi/3);
inter.coupling.matrix{1,1}=Q1;
inter.coupling.matrix{2,2}=Q2;

% Chemical shifts
inter.zeeman.scalar={0.0 3.0};

% Chemical exchange
inter.chem.parts={1,2};
inter.chem.rates=[-2e4   2e4
                   2e4  -2e4];
inter.chem.concs=[1.0 1.0];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'2H'};
parameters.rho0=state(spin_system,'L+','2H','chem');
parameters.coil=state(spin_system,'L+','2H');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=0.4e6;
parameters.npoints=1024;
parameters.zerofill=1024;
parameters.axis_units='Hz';
parameters.invert_axis=1;
parameters.grid='rep_2ang_200pts_sph';
parameters.rate=25000;
parameters.axis=[1 1 1];
parameters.max_rank=12;

% Acquisition
fid=singlerot(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'exp-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

