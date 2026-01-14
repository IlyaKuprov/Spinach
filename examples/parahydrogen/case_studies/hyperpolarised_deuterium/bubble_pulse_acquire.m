% Simulated PNL (partially negative line) spectrum of 
% ortho-deuterium in the presence of a parahydrogena-
% tion catalyst. Bubbling is followed by a 45-degree
% pulse. Paper link to follow in due course.
%
% thhu@mpinat.mpg.de
% anakin.aden@mpinat.mpg.de
% denismoll@hotmail.de
% julius.matz@mpinat.mpg.de
% stefan.gloeggler@mpinat.mpg.de
% ilya.kuprov@weizmann.ac.il

function bubble_pulse_acquire() 

% Spin system
sys.isotopes={'2H','2H','2H','2H'};

% Experimental chemical shifts
inter.zeeman.scalar={4.55 4.55 -13.5 -16.5};

% Experimental J-couplings
inter.coupling.scalar=cell(4,4);
inter.coupling.scalar{1,2} = 12.0;
inter.coupling.scalar{3,4} = 0.24;

% NQI tensors from a DFT calculation
inter.coupling.matrix{3,3}=1e3*[ 108.2    0.1   28.1
                                   0.1  -55.6    4.5
                                  28.1    4.5  -52.6];  
inter.coupling.matrix{4,4}=1e3*[ -55.2    8.1  -14.5
                                   8.1   -6.3  -73.9
                                 -14.5  -73.9   61.5];

% Cartesian coordinates, DFT calculation
inter.coordinates={[]; []; % None for D2
                   [-1.98  0.45  -0.55];
                   [-0.25  1.33  -1.65]};

% Kinetics
inter.chem.parts={[1 2],[3 4]};
inter.chem.rates=[-1  5000;...
                   1 -5000];
inter.chem.concs=[1 0];

% Magnet field
sys.magnet=7.05;

% Simulation formalsim
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'redfield','t1_t2'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={1e-12, 400e-12};
inter.r1_rates={0.04 0.04 0 0};
inter.r2_rates={8.00 8.00 0 0};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relaxation analysis
relaxan(spin_system);

% Spin states
[S,~,Q]=deut_pair(spin_system,1,2);

% Initial spin state: nothing
rho0=unit_state(spin_system);

% Assumptions
spin_system=assume(spin_system,'nmr');

% Evolution generators
H=hamiltonian(spin_system);
R=relaxation(spin_system);

% Free kinetics
KF=kinetics(spin_system);

% Kinetics with bubbling (a guess, needs proper rate)
pumped_state=S+Q{1}+Q{2}+Q{3}+Q{4}+Q{5}; pumped_state(1)=0;
KB=magpump(spin_system,KF,pumped_state,1e-1);

% Run the bubbling for 7 seconds
rho=evolution(spin_system,H+1i*R+1i*KB,[],rho0,7,1,'final');

% Pulse-acquire settings
parameters.spins={'2H'};
parameters.rho0=rho;
parameters.coil=state(spin_system,'L+',[1 2]);
parameters.pulse_op=operator(spin_system,'Ly','2H');
parameters.pulse_angle=pi/4;
parameters.decouple={};
parameters.offset=209.6554;
parameters.sweep=60;
parameters.npoints=256;
parameters.zerofill=1024;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@hp_acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); scale_figure([1.00 0.75]);
plot_1d(spin_system,real(spectrum),parameters);
kylabel('NMR intensity, a.u.'); xlim([4.40 4.70]);

end

