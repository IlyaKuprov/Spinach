% Evolution of state populations under ortho-deuterium
% bubbling in the presence of a parahydrogenation cata-
% lyst. No pulses, just evolution. Paper link to follow
% in due course.
%
% thhu@mpinat.mpg.de
% anakin.aden@mpinat.mpg.de
% denismoll@hotmail.de
% julius.matz@mpinat.mpg.de
% stefan.gloeggler@mpinat.mpg.de
% ilya.kuprov@weizmann.ac.il

function just_bubbling() 

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

% Pertinent spin states and coherences
[S,T,Q]=deut_pair(spin_system,1,2);

% Initial spin state: nothing
rho0=unit_state(spin_system);

% Detection: all relevant states, traceless
coils=[S    T{1} T{2} T{3} ...
       Q{1} Q{2} Q{3} Q{4} Q{5}]; 
coils(1,:)=0;

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

% Assemble the Liouvillians
LF=H+1i*R+1i*KF; % Free evolution
LB=H+1i*R+1i*KB; % Evolution with bubbling

% Run the bubbling for 7 seconds
traj_a=evolution(spin_system,LB,[],rho0,0.007,1000,'trajectory');

% Free evolution for 30 seconds
traj_b=evolution(spin_system,LF,[],traj_a(:,end),0.03,1000,'trajectory');

% Project out observables
obs=coils'*[traj_a traj_b];

% Plotting
time_axis=[linspace(0,7,1001) linspace(0,30,1001)+7];
kfigure(); scale_figure([1.00 0.75]);
plot(time_axis',real(obs([1 2 3 5 6 7],:)'));
kxlabel('time, seconds'); xlim tight; ylim padded;
kylabel('spin state coefficient, a.u.'); kgrid;
klegend('$S$','$T_{\pm1}$','$T_{0}$',...
        '$Q_{\pm2}$','$Q_{\pm1}$','$Q_{0}$');

end

