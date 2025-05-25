% A MAS DNP simulation performed as described in Fred Mentink-
% Vigier's paper (Spinach rotation conventions are different):
%
%         http://dx.doi.org/10.1016/j.jmr.2015.07.001
%
% Steady state rotor period simulation for a single crystal,
% computed using Newton-Raphson steady state solver.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function solid_effect_mas_steady()

% Magnet field
sys.magnet=9.403;

% Spin specification
sys.isotopes={'E','1H'};

% Interactions
inter.zeeman.eigs{1}=[2.00614  2.00194 2.00988];
inter.zeeman.euler{1}=pi*[253.6 105.1 123.8]/180;
inter.zeeman.eigs{2}=[0.00 0.00 0.00];
inter.zeeman.euler{2}=[0.00 0.00 0.00];
inter.coordinates={[0.00 0.00 0.00],...
                   [0.00 0.00 3.00]};

% Relaxation parameters
inter.relaxation={'weizmann'};
inter.weiz_r1e=1/0.3e-3;
inter.weiz_r1n=1/4.0;
inter.weiz_r2e=1/1.0e-6;
inter.weiz_r2n=1/0.2e-3;
inter.weiz_r1d=zeros(2,2);
inter.weiz_r2d=zeros(2,2);
inter.temperature=100;
inter.equilibrium='dibari';
inter.rlx_keep='secular';

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Stack generation parameters
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.rframes={};
parameters.orientation=[0 0 0];
parameters.spins={'E','1H'};
parameters.masframe='magnet';
parameters.offset=[0 0];
parameters.max_rank=3000;

% Rotor stack generation
H=rotor_stack(spin_system,parameters,'esr');

% Microwave operator
Hmw=operator(spin_system,'Lx','E');

% Zeeman offset operator
HzE=operator(spin_system,'Lz','E');

% Experiment parameters
parameters.mw_pwr=0.85e6;
parameters.mw_off=-400e6;
parameters.rate=12.5e3;

% Relaxation superoperator
R=relaxation(spin_system);

% Rotor period integration
nsteps=numel(H); P=speye(size(H{1}));
spin_system.sys.output='hush';
parfor n=1:nsteps %#ok<*PFBNS>
    P=propagator(spin_system,H{n}+2*pi*parameters.mw_pwr*Hmw+...
                                  2*pi*parameters.mw_off*HzE+1i*R,...
                                  1/(nsteps*parameters.rate))*P; 
end

% Steady state
rho_st=steady(spin_system,P,[],[],'newton');

% Rotor period trajectory
nsteps=numel(H); rho=zeros(numel(rho_st),nsteps,'like',1i); rho(:,1)=rho_st;
for n=2:nsteps
    rho(:,n)=step(spin_system,H{n}+2*pi*parameters.mw_pwr*Hmw+...
                                   2*pi*parameters.mw_off*HzE+1i*R,...
                                   rho(:,n-1),1/(nsteps*parameters.rate));
end

% Trajectory analysis
figure(); trajan(spin_system,rho,'level_populations');

% Thermal equilibrium state
rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'));

% Enhancement factor
Hz_dnp=mean(state(spin_system,'Lz','1H')'*rho);
Hz_eq=state(spin_system,'Lz','1H')'*rho_eq;
enh_factor=real(Hz_dnp/Hz_eq);
disp(['Enhancement factor: ' num2str(enh_factor)]);

end

