% Pulse-acquire FFT ESR version of the EasySpin Fremy salt test
% file, with acknowledgements to Stefan Stoll.
%
% The Spinach simulation is run using explicit time propagation
% in Liouville space with Redfield relaxation superoperator.
%
% Set to reproduce Figure 3a from
%
%           http://dx.doi.org/10.1209/epl/i2004-10459-y
%
% Calculation time: seconds
%
% matthew.krzystyniak@oerc.ox.ac.uk
% i.kuprov@soton.ac.uk

function relaxation_fremysalt()

% General layout
sys.magnet=0.33;
sys.isotopes={'E','14N'};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Interactions
inter.zeeman.eigs=cell(2,1);
inter.zeeman.euler=cell(2,1);
inter.zeeman.eigs{1}=[2.00785 2.00590 2.00265];
inter.zeeman.euler{1}=[0 0 0];
inter.coupling.eigs=cell(2,2);
inter.coupling.euler=cell(2,2);
inter.coupling.eigs{1,2}=[15.4137 14.0125 80.4316]*1e6;
inter.coupling.euler{1,2}=[0 0 0];

% Relaxation superoperator
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={8e-10};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E','cheap');
parameters.coil=state(spin_system,'L+','E','cheap');
parameters.decouple={};
parameters.offset=-2e7;
parameters.sweep=2e8;
parameters.npoints=512;
parameters.zerofill=1024;
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
xlabel('Electron Zeeman frequency, GHz');
ylabel('first derivative amplitude, a.u.');

end

