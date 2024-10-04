% 3D HCCH TOCSY experiment on a small protein fragment.
%
% Calculation time: seconds.
%
% m.walker@soton.ac.uk
% i.kuprov@soton.ac.uk

function hcch_tocsy_simple()

% Spin system
sys.isotopes={'13C','13C','13C','1H','1H'};
sys.labels={'C','CA','CB','HA','HB'};

% Interactions
sys.magnet=14.1;
inter.zeeman.scalar={180 60 40 4 2}; 
inter.coupling.scalar{2,3}=35.0;
inter.coupling.scalar{1,3}=0.3;
inter.coupling.scalar{1,2}=55.0;
inter.coupling.scalar{3,5}=130.0; 
inter.coupling.scalar{4,5}=7.0; 
inter.coupling.scalar{3,4}=6.0;
inter.coupling.scalar{2,4}=140.0;
inter.coupling.scalar{1,4}=4.0;
inter.coupling.scalar{5,5}=0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.J_ch=140;
parameters.delta=1.1e-3;
parameters.sl_tmix=2e-3;
parameters.lamp=1e4;
parameters.dipsi_dur=22.5e-3;
parameters.sweep=[4000 7500 4000];
parameters.spins={'1H','13C','1H'};
parameters.offset=[2000 7000 2000];
parameters.npoints=[64 64 64];
parameters.zerofill=[256 256 256];
parameters.decouple_f3={'13C'};
parameters.axis_units='ppm';

% Simulation
fid=liquid(spin_system,@hcch_tocsy,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'cos'},{'cos'},{'cos'}});

% Fourier transform
spectrum=fftshift(fftn(fid,parameters.zerofill));

% Plotting
figure(); plot_3d(spin_system,abs(spectrum),parameters,...
                  10,[0.05 0.25 0.05 0.25],2,'positive');

end

