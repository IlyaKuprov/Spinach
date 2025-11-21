% A model of the CIDNP magnetisation pumping process described in
% IK's paper:
% 
%           https://doi.org/10.1016/j.jmr.2004.01.011
%
% The system is pumped for 0.5 seconds, and then allowed to relax.
%
% Calculation time: seconds.
%
% Miguel Mompean
% Ilya Kuprov

function cidnp_flash_acquire()

% Magnet field
sys.magnet=14.1;

% Isotopes
sys.isotopes={'1H','19F'};

% Chemical shifts
inter.zeeman.scalar={0.0, 0.0};

% Chemical shift anisotropies (DFT)
inter.zeeman.eigs{1}=[0 0 0];
inter.zeeman.euler{1}=[0 0 0];
inter.zeeman.eigs{2}=[-47 -16  63]; 
inter.zeeman.euler{2}=[0 0 0];

% Coordinates (DFT)
inter.coordinates={[0.00 0.00 0.00],...
                   [0.00 2.60 0.00]};

% J-coupling (expt)
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=50;

% Relaxation theories
inter.relaxation={'redfield'};
inter.rlx_keep='secular';
inter.tau_c={110e-12};
inter.equilibrium='zero';

% Formalism and basis
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get the relevant states
Hz=state(spin_system,{'Lz'},{1});
Fz=state(spin_system,{'Lz'},{2});
HzFz=state(spin_system,{'Lz','Lz'},{1,2});
Hz=Hz/norm(Hz,2); Fz=Fz/norm(Fz,2);
HzFz=HzFz/norm(HzFz,2); 

% Get the Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Get relaxation matrix
R=relaxation(spin_system);

% Thermalise R to Hz+Fz
R(1,1)=1; R(:,1)=-R*(Hz+Fz);

% Add pumping terms
R_light=magpump(spin_system,R,Hz,1.3);
R_light=magpump(spin_system,R_light,Fz,34.0);

% Start at equilibrium
rho0=unit_state(spin_system)+Hz+Fz;

% Run the evolution with illumination, 50 steps of 0.01 seconds
traj_a=evolution(spin_system,H+1i*R_light,[],rho0,0.01,50,'trajectory');

% Run the evolution with illumination, 500 steps of 0.01 seconds
traj_b=evolution(spin_system,H+1i*R,[],traj_a(:,end),0.01,500,'trajectory'); 

% Concatenate trajectories
traj=[traj_a traj_b(:,2:end)];
                          
% Get the time axis
x_axis=linspace(0,5.5,551);

% Do the plotting
kfigure(); scale_figure([1.5 1.0]);
subplot(1,3,1); plot(x_axis,Fz'*traj); kgrid;
ktitle('$F_{\mathrm{Z}}$'); kxlabel('time, s'); 
axis tight; ylim([0 15]);
subplot(1,3,2); plot(x_axis,Hz'*traj); kgrid;
ktitle('$H_{\mathrm{Z}}$'); kxlabel('time, s'); 
axis tight; ylim([0 1.5]);
subplot(1,3,3); plot(x_axis,-2*HzFz'*traj); kgrid;
ktitle('$H_{\mathrm{Z}}F_{\mathrm{Z}}$'); 
kxlabel('time, s'); axis tight; ylim([-1.0 0]);

end

