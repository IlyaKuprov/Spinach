% A minimal example of HNCA pulse sequence simulation.
%
% Calculation time: seconds.
%
% m.walker@soton.ac.uk
% i.kuprov@soton.ac.uk

function hnca_simple()

% Magnet field
sys.magnet=14.1;

% Spin system
sys.isotopes={'15N','13C','1H','13C'};
sys.labels={'N','CA','H','C'};

% Interactions
inter.zeeman.scalar={110 60 8 180};
inter.coupling.scalar=cell(4);
inter.coupling.scalar{1,3}=92;
inter.coupling.scalar{1,2}=11;
inter.coupling.scalar{1,4}=15;
inter.coupling.scalar{2,4}=55;
inter.coupling.scalar{2,3}=2;
inter.coupling.scalar{3,4}=4;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H','13C','15N'};
parameters.sweep=[3000 5000 2800];
parameters.offset=[5100 8600 -7200];
parameters.npoints=[64 64 64];
parameters.zerofill=[256 256 256];
parameters.axis_units='ppm';

% Simulation
fid=liquid(spin_system,@hnca,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'sqcos'},{'sqcos'},{'sqcos'}});

% Fourier transform
spectrum=fftshift(fftn(fid,parameters.zerofill));

% Plotting
figure(); plot_3d(spin_system,abs(spectrum),parameters,...
                  10,[0.2 0.9 0.2 0.9],2,'positive');

end

