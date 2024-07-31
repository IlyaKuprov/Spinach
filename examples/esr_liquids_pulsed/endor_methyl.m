% Mims ENDOR spectrum of a methyl radical in liquid state. 
% Magnetic parameters taken from a DFT calculation.
%
% Calculation time: seconds
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk

function endor_methyl()

% Ignore coordinate information (HFCs provided)
options.no_xyz=1;

% System properties (vacuum DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/methyl.log'),...
                         {{'E','E'},{'H','1H'}},[0 0],options);
% Magnet field
sys.magnet=0.33;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.sym_group={'S3'};
bas.sym_spins={[1 2 3]};

% Sequence parameters
parameters.offset=0;
parameters.npoints=512;
parameters.sweep=4e8;
parameters.tau=100e-9;
parameters.zerofill=4096;
parameters.spins={'E'};
parameters.axis_units='MHz';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@endor_mims,parameters,'esr');

% Crude apodization
fid=apodization(fid-mean(fid),'kaiser-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,abs(spectrum),parameters);

end

