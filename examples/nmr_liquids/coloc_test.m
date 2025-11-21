% A simple COLOC pulse sequence example for a two-spin
% 1H-13C system with a long-range J-coupling.
%
% Calculation time: seconds.
%
% b.macaulay@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function coloc_test()

% Magnet field
sys.magnet=11.7;

% Spin system
sys.isotopes={'1H','13C'};

% Interactions
inter.zeeman.scalar={4.0 75.0};
inter.coupling.scalar{1,2}=5.0;
inter.coupling.scalar{2,2}=0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Sequence parameters
parameters.spins={'1H','13C'};
parameters.delta2=30e-3;
parameters.offset=[2250 5000];
parameters.sweep=[5000 12000];
parameters.npoints=[256 256];
parameters.zerofill=[512 512];
parameters.axis_units='ppm';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@coloc,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'cos'},{'cos'}});
    
% Fourier transform
spec=fftshift(fft2(fid,parameters.zerofill(2),...
                       parameters.zerofill(1)));
% Plot
kfigure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,abs(spec),parameters,...
        20,[0.1 0.5 0.1 0.5],2,256,6,'both');

end

