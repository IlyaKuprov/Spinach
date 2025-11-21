% CT-COSY of three spin system.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il
% mrw1g16@soton.ac.uk

function ct_cosy_three_spin()

% Spin system and interactions
sys.magnet=14.1;
sys.isotopes={'1H','1H','1H'};
inter.zeeman.scalar={2.70 4.10 6.50};
inter.coupling.scalar{1,2}=10;
inter.coupling.scalar{2,3}=8;
inter.coupling.scalar{1,3}=4;
inter.coupling.scalar{3,3}=0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=4.0;

% Sequence parameters
parameters.offset=2700;
parameters.sweep=[3500 3500];
parameters.npoints=[256 256];
parameters.zerofill=[512 512];
parameters.spins={'1H'};
parameters.axis_units='ppm';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@ct_cosy,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'sqcos'},{'sqcos'}}); 

% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill(2),...
                           parameters.zerofill(1)));

% Plotting
kfigure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,abs(spectrum),parameters,...
        20,[0.05 0.25 0.05 0.25],2,256,6,'positive');

end

