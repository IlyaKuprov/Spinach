% A time-domain pulse-acquire version of the EasySpin biaryl test file, 
% with acknowledgements to Stefan Stoll. 
% 
% The Spinach simulation is run using explicit time propagation in Liou-
% ville space. Full symmetry treatment using the S2xS2xS2xS2xS2xS2 group
% direct product is performed.
%
% Calculation time: seconds
%
% matthew.krzystyniak@oerc.ox.ac.uk
% ilya.kuprov@weizmann.ac.il

function pulse_acquire_biaryl()

% Magnet induction
sys.magnet=0.33;

% Isotope list
sys.isotopes={'E','14N','1H','1H','1H','1H','1H',...
                  '14N','1H','1H','1H','1H','1H'};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.longitudinals={'1H','14N'};
bas.projections=+1;
bas.sym_group={'S2','S2','S2','S2','S2','S2'};
bas.sym_spins={[2 8],[3 9],[4 10],[5 11],[6 12],[7 13]};

% Zeeman interactions
inter.zeeman.scalar=cell(13,1);
inter.zeeman.scalar{1}=2.00316;

% Spin-spin couplings
inter.coupling.scalar=cell(13,13);
inter.coupling.scalar{1,2} =  12.16e6;
inter.coupling.scalar{1,3} =  -6.7e6;
inter.coupling.scalar{1,4} =  -1.82e6;
inter.coupling.scalar{1,5} =  -7.88e6;
inter.coupling.scalar{1,6} =  -0.64e6;
inter.coupling.scalar{1,7} =  67.93e6;
inter.coupling.scalar{1,8} =  12.16e6;
inter.coupling.scalar{1,9} =  -6.7e6;
inter.coupling.scalar{1,10} = -1.82e6;
inter.coupling.scalar{1,11} = -7.88e6;
inter.coupling.scalar{1,12} = -0.64e6;
inter.coupling.scalar{1,13} = 67.93e6;

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=5e5;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E');
parameters.coil=state(spin_system,'L+','E');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=3e8;
parameters.npoints=4096;
parameters.zerofill=16384;
parameters.axis_units='GHz-labframe';
parameters.derivative=1;
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'esr');

% Apodisation
fid=apodisation(spin_system,fid,{{'none'}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plot the spectrum
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

