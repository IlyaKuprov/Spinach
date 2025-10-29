% Bright fat effect under CPMG echo train - magnetisation losses
% are greater in MRI experiments on J-coupled systems because co-
% herences are lost in the depths of the Hilbert space.
%
% Simulation time: minutes, faster with a Tesla V100 GPU.
%
% a.j.allami@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function bright_fat_effect_cpmg()

% Magnetic induction
sys.magnet=3.0;

% Spin system
sys.isotopes={'1H','1H','1H','1H','1H','1H'};

% Chemical shifts
inter.zeeman.scalar={1.0, 2.0, 3.0, 1.0, 2.0, 3.0};

% J-coupling
inter.coupling.scalar=cell(6,6);
inter.coupling.scalar{1,2}=0;
inter.coupling.scalar{2,3}=0;
inter.coupling.scalar{1,3}=0;
inter.coupling.scalar{4,5}=11;
inter.coupling.scalar{4,6}=17;
inter.coupling.scalar{5,6}=23;

% Spins 1,2,3 are molecule A
% and 4,5,6 are molecule B
inter.chem.parts={[1 2 3],[4 5 6]};

% Kinetic rate matrix (Hz)
inter.chem.rates=[0 0; 0 0];
inter.chem.concs=[1 1];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Disable path tracing
sys.disable={'pt'};

% This needs a GPU
% sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.decouple={};
parameters.offset=0.0;
parameters.npulses=48;
parameters.dec_time=80e-3;

% Sample geometry
parameters.dims=[0.30 0.25];
parameters.npts=[100 200];
parameters.deriv={'period',3};

% Relaxation phantoms and operators
parameters.rlx_ph={}; parameters.rlx_op={};

% Initial and detection state phantoms
load('../../etc/phantoms/bright_fat_left.mat','left');  
load('../../etc/phantoms/bright_fat_right.mat','right'); 
parameters.rho0_ph={1-left,1-right};
parameters.rho0_st={state(spin_system,'Lz',[1 2 3]),...
                    state(spin_system,'Lz',[4 5 6])};
parameters.coil_ph={ones(parameters.npts)};
parameters.coil_st={state(spin_system,'Lx','1H')};

% No diffusion or flow
parameters.u=zeros(parameters.npts);
parameters.v=zeros(parameters.npts);
parameters.diff=0;

% Run the simulation
mri=imaging(spin_system,@cpmg_dec,parameters);

% Plotting
figure(); surf(abs(mri)); set(gca,'XDir','reverse');
kxlabel('FOV1, px'); kylabel('FOV2, px'); kgrid;
ktitle('Bright fat effect under CPMG echo train');

end

