% A pulse-acquire FFT version of the EasySpin parafluoronitrobenzene test
% file, with acknowledgements to Stefan Stoll.
%
% The Spinach simulation is run using explicit time propagation in Liouville
% space, including secular Redfield relaxation superoperator.
%
% Calculation time: minutes
%
% matthew.krzystyniak@oerc.ox.ac.uk
% i.kuprov@soton.ac.uk

function relaxation_parafluoronitrobenzene()

% Magnet field
sys.magnet=0.33898;

% Isotope list
sys.isotopes={'E','14N','19F','1H','1H','1H','1H'};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.zero_quantum={'1H'};
bas.sym_group={'S2','S2'};
bas.sym_spins={[4 5],[6 7]};

% Zeeman interactions
inter.zeeman.eigs=cell(7,1);
inter.zeeman.euler=cell(7,1);
inter.zeeman.eigs{1}=[2.0032 2.0012 2.0097];
inter.zeeman.euler{1}=[0 0 0];

% Spin-spin couplings
inter.coupling.eigs=cell(7,7);
inter.coupling.euler=cell(7,7);
inter.coupling.eigs{1,2}=(40.40+[24 -12 -12])*1e6;
inter.coupling.eigs{1,3}=(22.51+[34.9 -19.8 -15])*1e6;
inter.coupling.eigs{1,4}=[9.69 9.69 9.69]*1e6;
inter.coupling.eigs{1,5}=[9.69 9.69 9.69]*1e6;
inter.coupling.eigs{1,6}=[3.16 3.16 3.16]*1e6;
inter.coupling.eigs{1,7}=[3.16 3.16 3.16]*1e6;
inter.coupling.euler{1,2}=[0 0 0];
inter.coupling.euler{1,3}=[0 0 0];
inter.coupling.euler{1,4}=[0 0 0];
inter.coupling.euler{1,5}=[0 0 0];
inter.coupling.euler{1,6}=[0 0 0];
inter.coupling.euler{1,7}=[0 0 0];

% Relaxation superoperator
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={160e-12};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E','cheap');
parameters.coil=state(spin_system,'L+','E','cheap');
parameters.decouple={};
parameters.offset=-1e7;
parameters.sweep=2e8;
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

