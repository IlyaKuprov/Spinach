% Figure 8.26 from Andrew Derome's "Modern NMR Techniques
% for Chemistry Research".
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il
% gareth.charnock@oerc.ox.ac.uk

function cosy90_derome()

% Magnet field
sys.magnet=16.1;

% Spin system and interactions
sys.isotopes={'1H','1H','1H'};
inter.zeeman.scalar={3.70 3.92 4.50};
inter.coupling.scalar{1,2}=10;
inter.coupling.scalar{2,3}=12;
inter.coupling.scalar{1,3}=4;
inter.coupling.scalar{3,3}=0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=4.0;

% Sequence parameters
parameters.angle=pi/2;
parameters.offset=2800;
parameters.sweep=700;
parameters.npoints=[1024 1024];
parameters.zerofill=[2048 2048];
parameters.spins={'1H'};
parameters.axis_units='ppm';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@cosy,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'sqcos'},{'sqcos'}});

% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill(2),...
                           parameters.zerofill(1)));

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,real(spectrum),parameters,...
        20,[0.1 0.5 0.1 0.5],2,256,6,'both');

end

