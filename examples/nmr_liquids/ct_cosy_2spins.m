% CT COSY spectrum for 2 spins.
% 
% http://dx.doi.org/10.1002/jhet.5570250160
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% mrw1g16@soton.ac.uk

function ct_cosy_2spins()

% Spin system
sys.isotopes={'1H','1H'};

% Interactions
sys.magnet=5.9;
inter.zeeman.scalar={1.00 3.00};
inter.coupling.scalar{1,2}=7.0; 
inter.coupling.scalar{2,2}=0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Sequence parameters
parameters.offset=500;
parameters.sweep=[2000 2000];
parameters.npoints=[512 512];
parameters.zerofill=[2048 2048];
parameters.spins={'1H'};
parameters.axis_units='ppm';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@ct_cosy,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'cos'},{'cos'}}); 

% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill(2),...
                           parameters.zerofill(1)));

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,abs(spectrum),parameters,...
        20,[0.1 0.5 0.1 0.5],2,256,6,'positive');

end

