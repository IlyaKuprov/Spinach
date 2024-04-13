% Powder-averaged HYSCORE on a 14N nitroxide radical. Time-domain
% simulation in Liouville space. Set to reproduce Figure 2a from 
%
%           http://dx.doi.org/10.1080/00268979809483260
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function hyscore_nitroxide_powder()

% Magnet field
sys.magnet=0.350;

% System specification
sys.isotopes={'14N','E'};
inter.zeeman.scalar={0.0 2.0000};
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=5e6;
inter.coupling.matrix=cell(2,2);
inter.coupling.matrix{1,1}=eeqq2nqi(2.4e6,0.5,1,[0 0 0]);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Disable trajectory-level SSR algorithms
sys.disable={'trajlevel','colorbar'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set the sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'Lz','E');
parameters.coil=state(spin_system,'L+','E');
parameters.offset=0;
parameters.nsteps=[128 128];
parameters.zerofill=[256 256];
parameters.sweep=20e6;
parameters.tau=136e-9;
parameters.grid='rep_2ang_800pts_sph';
parameters.axis_units='MHz';
parameters.verbose=0;

% Simulation
fid=powder(spin_system,@hyscore,parameters,'esr');

% Centre signal suppression
fid=fid-mean(mean(fid));

% Apodization
fid=apodization(fid,'cosbell-2d');

% Fourier transform
spectrum=fftshift(fft2(fid,parameters.zerofill(2),parameters.zerofill(1)));

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,abs(spectrum),parameters,20,...
        [0.025 0.5 0.025 0.5],2,256,6,'positive');

end

