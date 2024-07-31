% CW ENDOR on 2-methoxy-1,4-benzoquinone radical in liquid state. Set to
% reproduce Figure 2 in http://dx.doi.org/10.1002/mrc.1260280313
%
% Calculation time: seconds
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk

function endor_benzoquinone()

% Magnet field
sys.magnet=0.33;

% Isotopes and interactions
sys.isotopes={'E','1H','1H','1H','1H','1H','1H'};
inter.zeeman.scalar={2.004577 0 0 0 0 0 0};
inter.coupling.scalar=cell(7,7);
inter.coupling.scalar{2,1}=mt2hz(0.08); 
inter.coupling.scalar{3,1}=mt2hz(0.08); 
inter.coupling.scalar{4,1}=mt2hz(0.08); 
inter.coupling.scalar{5,1}=mt2hz(-0.059); 
inter.coupling.scalar{6,1}=mt2hz(-0.364); 
inter.coupling.scalar{7,1}=mt2hz(-0.204); 

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.sym_group={'S3'};
bas.sym_spins={[2 3 4]};

% Sequence parameters
parameters.offset=0;
parameters.sweep=50e6;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.spins={'E'};
parameters.axis_units='MHz';
parameters.derivative=1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@endor_cw,parameters,'esr');

% Crude apodization
fid=apodization(fid-mean(fid),'kaiser-1d',20);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,-abs(spectrum),parameters);

end

