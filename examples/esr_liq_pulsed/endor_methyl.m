% Mims ENDOR spectrum of a methyl radical in liquid state. 
% Magnetic parameters taken from a DFT calculation.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function endor_methyl()

% Magnet field
sys.magnet=0.346;

% Spin system and interactions
sys.isotopes={'1H','1H','1H','E'};
inter.zeeman.scalar={0 0 0 2.0026};
inter.coupling.scalar=cell(4,4);
inter.coupling.scalar{1,4}=-64.4e6;
inter.coupling.scalar{2,4}=-64.4e6;
inter.coupling.scalar{3,4}=-64.4e6;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.sym_group={'S3'};
bas.sym_spins={[1 2 3]};

% Sequence parameters
parameters.offset=0;
parameters.npoints=512;
parameters.sweep=1.2e8;
parameters.tau=100e-9;
parameters.zerofill=4096;
parameters.spins={'E'};
parameters.axis_units='MHz';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@endor_mims,parameters,'esr');

% Crude apodisation
fid=apodisation(spin_system,fid-mean(fid),{{'kaiser',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,abs(spectrum),parameters);
kxlabel('Nuclear frequency, MHz');

end

