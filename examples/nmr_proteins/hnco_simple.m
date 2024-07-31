% A minimal example of HNCO pulse sequence simulation.
%
% Calculation time: seconds.
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk

function hnco_simple()

% Magnet field
sys.magnet=14.1;

% Spin system
sys.isotopes={'15N','13C','1H','13C'};
sys.labels={'N','CA','H','C'};

% Interactions
inter.zeeman.scalar={110 55 8.0 180};
inter.coupling.scalar=cell(4);
inter.coupling.scalar{1,3}=92;
inter.coupling.scalar{1,2}=11;
inter.coupling.scalar{1,4}=15;
inter.coupling.scalar{2,4}=55;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'15N','13C','1H'};
parameters.offset=[-7200 25000 4800];
parameters.sweep=[5000 10000 5000];
parameters.npoints=[63 64 65];
parameters.zerofill=[255 256 257];
parameters.tau=[2.25e-3 14e-3 4e-3];
parameters.f1_decouple=1;
parameters.axis_units='ppm';

% Simulation
fid=liquid(spin_system,@hnco,parameters,'nmr');

% Apodization
fid=apodization(fid,'sqcosbell-3d');

% Fourier transform
spectrum=fftn(fid,parameters.zerofill);
spectrum=fftshift(spectrum);

% Plotting
figure(); plot_3d(spin_system,abs(spectrum),parameters,...
                  10,[0.2 0.9 0.2 0.9],2,'positive');

end

