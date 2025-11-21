% Magnitude mode COSY-45 spectrum of rotenone.
% 
% http://dx.doi.org/10.1002/jhet.5570250160
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% gareth.charnock@oerc.ox.ac.uk

function cosy45_rotenone()

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H',...
              '1H','1H','1H','1H','1H','1H','1H','1H','1H',...
              '1H','1H','1H','1H'};

% Interactions
sys.magnet=5.9;
inter.zeeman.scalar={6.72 6.40 4.13 4.56 4.89 6.46 7.79 3.79 2.91...
                     3.27 5.19 4.89 5.03 1.72 1.72 1.72 3.72 3.72...
                     3.72 3.76 3.76 3.76};
inter.coupling.scalar{3,4}=12.1; 
inter.coupling.scalar{4,5}=3.1; 
inter.coupling.scalar{3,5}=1.0; 
inter.coupling.scalar{3,8}=1.0; 
inter.coupling.scalar{1,8}=1.0;
inter.coupling.scalar{6,7}=8.6; 
inter.coupling.scalar{5,8}=4.1; 
inter.coupling.scalar{7,9}=0.7; 
inter.coupling.scalar{7,10}=0.7; 
inter.coupling.scalar{9,10}=15.8;
inter.coupling.scalar{10,11}=9.8; 
inter.coupling.scalar{9,11}=8.1; 
inter.coupling.scalar{13,14}=1.5; 
inter.coupling.scalar{12,14}=0.9; 
inter.coupling.scalar{22,22}=0;

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.space_level=1;
bas.connectivity='scalar_couplings';
bas.sym_group={'S3','S3','S3'};
bas.sym_spins={[14 15 16],[17 18 19],[20 21 22]};

% Sequence parameters
parameters.angle=pi/4;
parameters.offset=1200;
parameters.sweep=2000;
parameters.npoints=[512 512];
parameters.zerofill=[2048 2048];
parameters.spins={'1H'};
parameters.axis_units='ppm';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@cosy,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'cos'},{'cos'}});

% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill(2),...
                           parameters.zerofill(1)));

% Plotting
kfigure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,abs(spectrum),parameters,...
        20,[0.0025 0.05 0.0025 0.05],2,256,6,'positive');

end

