% Pulse-acquire FFT ESR on 2-methoxy-1,4-benzoquinone radical in 
% liquid state. Set to reproduce Figure 1 in
%
%             http://dx.doi.org/10.1002/mrc.1260280313
%
% Simple common linewidth is used as a relaxation model.
%
% Calculation time: seconds
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk

function pulse_acquire_benzoquinone()

% Magnet induction
sys.magnet=0.33;

% Isotope list
sys.isotopes={'E','1H','1H','1H','1H','1H','1H'};

% Zeeman interactions and couplings
inter.zeeman.scalar={2.004577 0 0 0 0 0 0};
inter.coupling.scalar=cell(7,7);
inter.coupling.scalar{2,1}=mt2hz(0.08); 
inter.coupling.scalar{3,1}=mt2hz(0.08); 
inter.coupling.scalar{4,1}=mt2hz(0.08); 
inter.coupling.scalar{5,1}=mt2hz(-0.059); 
inter.coupling.scalar{6,1}=mt2hz(-0.364); 
inter.coupling.scalar{7,1}=mt2hz(-0.204); 

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=1e6;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.longitudinals={'1H'};
bas.projections=+1;
bas.sym_group={'S3'};
bas.sym_spins={[2 3 4]};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E','cheap');
parameters.coil=state(spin_system,'L+','E','cheap');
parameters.decouple={};
parameters.offset=-1e7;
parameters.sweep=3e7;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.axis_units='GHz-labframe';
parameters.derivative=1;
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'esr');

% Apodization
fid=apodization(fid,'none-1d');

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

