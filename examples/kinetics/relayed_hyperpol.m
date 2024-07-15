% Relayed NOE from hyperpolarized water to ALA-GLY dipeptide,
% generating Figure S7 from 
%
%          https://doi.org/10.1016/j.jmr.2024.107727
% 
% Christopher Pötzl

function relayed_hyperpol()

% Simulation timing parameters
dt=0.125; nsteps=128;

% Magnet field
sys.magnet=16.4;

% 30 protons in the system
sys.isotopes=repelem({'1H'},30);

% Cartesian coordinates of pertinent protons
inter.coordinates={[6.67  4.45  4.03];  % labile
                   [6.40  3.10  4.83];  % labile
                   [7.42  4.39  5.52];  % labile
                   [4.64  5.97  3.40];  % labile
                   [4.50  4.35  5.28];  % aliphatic 
                   [6.44  5.29  7.63];  % aliphatic
                   [5.77  3.57  7.58];  % aliphatic
                   [4.70  4.88  7.79];  % aliphatic
                   [5.33  8.70  4.21];  % aliphatic
                   [4.44  8.22  2.75]}; % aliphatic

% 20 water protons exist, but have no coordinates to
% prevent direct cross-relaxation from happening
inter.coordinates=[inter.coordinates; repelem({[]},20)'];

% Chemical shifts, all water at 4.5 ppm
inter.zeeman.scalar={8.45 8.45 8.45 8.11 3.73 ...
                     0.99 0.99 0.99 3.99 3.32};
inter.zeeman.scalar=[inter.zeeman.scalar repelem({4.5},20)];

% Relaxation theories
inter.relaxation={'redfield','t1_t2'};
inter.equilibrium='dibari';
inter.rlx_keep='secular';
inter.tau_c={1.2e-10};
inter.temperature=298;

% Empirical relaxation at 0.1 Hz for water
inter.r1_rates=num2cell([zeros(1,10) 0.1*ones(1,20)]);
inter.r2_rates=num2cell([zeros(1,10) 0.1*ones(1,20)]);

% Basis set, single-spin for water, up to 
% three-spin orders for the molecule
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='full_tensors';
bas.space_level=3;
bas.level=1;

% Exchange flux matrix
inter.chem.flux_rate=zeros(30,30);
inter.chem.flux_rate(1:4,11:20)=20;
inter.chem.flux_rate(11:20,1:4)=20;
inter.chem.flux_type='intermolecular';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system=assume(spin_system,'nmr');
        
% Build the Liouvillian
H=hamiltonian(spin_system);
R=relaxation(spin_system);
K=kinetics(spin_system);
L=H+1i*R+1i*K;
       
% Get equilibrium density matrix
HL=hamiltonian(assume(spin_system,'labframe'),'left');
rho=equilibrium(spin_system,HL);

% Polarise the water 100%
Wz=state(spin_system,'Lz',11:20);
rho=rho-Wz*(Wz'*rho)/norm(Wz,2)^2+Wz;

% Get detection states
H_aliph=[6 7 8]; H_alpha=5;
HZ_aliph=state(spin_system,'Lz',H_aliph); 
HZ_alpha=state(spin_system,'Lz',H_alpha);
            
% Time evolution simulation
result=evolution(spin_system,L,[HZ_aliph HZ_alpha],...
                 rho,dt,nsteps,'multichannel');
    
% Plotting
figure('Name','HA and CH3 magnetization evolution');
time_axis=linspace(0,nsteps*dt,nsteps+1); 
plot(time_axis,real(result)); xlim tight; kgrid;
kxlabel('time, seconds'); kylabel('magnetisation, a.u.'); 
klegend({'$CH_{3}$','$H_{\alpha}$'},'Location','SouthEast'); 

end

