% Relayed NOE from hyperpolarized water to ALA-GLY dipeptide,
% generating Figure S7 from 
%
%          https://doi.org/10.1016/j.jmr.2024.107727
% 
% Additive replacement records retain internal peptide orders and discard
% cross-molecule orders, as in the intermolecular exchange model. Every
% pool has invariant unit concentration, permitting a constant generator.
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

% Coordinate-free water has no direct cross-relaxation
inter.coordinates=[inter.coordinates; repelem({[]},20)'];

% Chemical shifts, all water at 4.5 ppm
inter.zeeman.scalar={8.45 8.45 8.45 8.11 3.73 ...
                     0.99 0.99 0.99 3.99 3.32};
inter.zeeman.scalar=[inter.zeeman.scalar repelem({4.5},20)];

% Relaxation theories
inter.relaxation={'redfield','t1_t2'};
inter.equilibrium='dibari';
inter.rlx_keep='secular';
inter.tau_c=repmat({1.2e-10},1,21);
inter.temperature=298;

% Empirical relaxation at 0.1 Hz for water
inter.r1_rates=num2cell([zeros(1,10) 0.1*ones(1,20)]);
inter.r2_rates=num2cell([zeros(1,10) 0.1*ones(1,20)]);

% Peptide and twenty independent unit-concentration water pools
inter.chem.parts=[{1:10} num2cell(11:30)];
inter.chem.concs=ones(1,21);

% Intermolecular replacements preserve populations and trace departing spins
inter.chem.reactions=cell(1,40);
for n=1:4
    for k=11:20
        matching=[(1:10)' (1:10)';k k];
        matching([n 11],2)=[k;n];
        inter.chem.reactions{10*(n-1)+k-10}=...
            struct('reactants',[1 k-9],'products',[1 k-9],...
                   'matching',matching,'rate',20,'closure','additive');
    end
end

% Three-spin peptide orders and complete one-spin water bases
bas.formalism='sphten-liouv';
bas.approximation=repmat({'IK-1'},1,21);
bas.connectivity=repmat({'full_tensors'},1,21);
bas.prox_level=[{3} repmat({1},1,20)];
bas.inter_level=repmat({1},1,21);

% Enable zero track elimination
sys.enable={'zte'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system=assume(spin_system,'nmr');
        
% Freeze additive chemistry at its invariant unit concentrations
H=hamiltonian(spin_system);
R=relaxation(spin_system);
K=kinetics(spin_system);
K=K(0,unit_state(spin_system));
L=H+1i*R+1i*K;
       
% Isotropic thermal equilibrium
rho=equilibrium(spin_system);

% Polarise the water 100%
Wz=coil_state(spin_system,'Lz',11:20,'exact');
rho=rho-Wz*(Wz'*rho)/norm(Wz,2)^2+Wz;

% Get detection states
H_aliph=[6 7 8]; H_alpha=5;
HZ_aliph=coil_state(spin_system,'Lz',H_aliph,'exact');
HZ_alpha=coil_state(spin_system,'Lz',H_alpha,'exact');
            
% Time evolution simulation
result=evolution(spin_system,L,[HZ_aliph HZ_alpha],...
                 rho,dt,nsteps,'multichannel');
    
% Plotting
time_axis=linspace(0,nsteps*dt,nsteps+1); 
figure('Name','HA and CH3 magnetisation evolution');
scale_figure([1.00 0.65]); plot(time_axis,real(result)); 
kxlabel('time, seconds'); xlim tight; kgrid;
kylabel('magnetisation, a.u.'); ylim padded; 
klegend({'CH$_{3}$','H$_{\alpha}$'},'Location','Best');

end

