% Lucio Frydman's water exchange based spin-lock pump, Figure 9
% from https://doi.org/10.1016/j.jmr.2021.107083
%
% Calculation time: seconds
%
% mihajlo.novakovic@weizmann.ac.il
% lucio.frydman@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il
% mariagrazia.concilio@sjtu.edu.cn

function frydman_pump_b()

% Number of water protons
n_water_protons=100;

% Magnet field 
sys.magnet=11.7;

% Core spin system 
sys.isotopes={'1H','15N','13C','13C'};

% Add water protons
sys.isotopes(5:(5+n_water_protons-1))={'1H'};

% Chemical shifts
inter.zeeman.scalar={0.0,0.0,0.0,0.0};
inter.zeeman.scalar(5:(5+n_water_protons-1))={0.0};

% Scalar couplings
inter.coupling.scalar=cell(n_water_protons+4);
inter.coupling.scalar{1,2}=-45;  
inter.coupling.scalar{2,4}=8; 

% Estimated relaxation times, seconds
T1H=0.2722; T1N=0.8; T1Ca=2; T1CO=2; T1W=0.2994;
T2H=0.2722; T2N=0.8; T2Ca=2; T2CO=2; T2W=0.2994;

% Relaxation theory
inter.relaxation={'t1_t2'};
inter.r1_rates=num2cell(1./[T1H T1N T1Ca T1CO T1W*ones(1,n_water_protons)]);
inter.r2_rates=num2cell(1./[T2H T2N T2Ca T2CO T2W*ones(1,n_water_protons)]);
inter.rlx_keep='diagonal';
inter.equilibrium='IME';
inter.temperature=298;

% Peptide and independent water-pool substances
inter.chem.parts=[{1:4} num2cell(5:(n_water_protons+4))];
inter.chem.concs=ones(1,n_water_protons+1);

% Per-substance basis settings retain each internal spin space
bas.formalism='sphten-liouv';
bas.approximation=repmat({'IK-1'},1,n_water_protons+1);
bas.connectivity=repmat({'scalar_couplings'},1,n_water_protons+1);
bas.inter_level=[{4} repmat({1},1,n_water_protons)];
bas.prox_level=repmat({1},1,n_water_protons+1);

% Exchange rates between NH and nearest water, and among water protons
nh_wt_exch_rate=1000; wt_wt_exch_rate=1e4;

% Replacement swaps the amide proton while retaining the peptide core
inter.chem.reactions={struct('reactants',[1 2],'products',[1 2],...
                            'matching',[1 5;2 2;3 3;4 4;5 1],...
                            'rate',nh_wt_exch_rate)};

% Every water pair exchanges once as a symmetric replacement event
for n=1:n_water_protons
    for k=(n+1):n_water_protons
        inter.chem.reactions{end+1}=struct('reactants',[n+1 k+1],...
            'products',[n+1 k+1],'matching',[n+4 k+4;k+4 n+4],...
            'rate',wt_wt_exch_rate);
    end
end

% Enable zero track elimination
sys.enable={'zte'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation parameters
parameters.spins={'1H','15N','13C'};
parameters.cp_dur=[11e-3 53e-3];
parameters.cp_npt=[11 53];
parameters.nloops=10;

% Get system trajectory (sequence is below)
traj=liquid(spin_system,@frydman_pump,parameters,'nmr');

% Observables: peptide bond spins
Hz=coil_state(spin_system,{'Lz'},{1},'exact');
Hx=coil_state(spin_system,{'Lx'},{1},'exact');
Nz=coil_state(spin_system,{'Lz'},{2},'exact');
Nx=coil_state(spin_system,{'Lx'},{2},'exact');
COz=coil_state(spin_system,{'Lz'},{4},'exact');
COx=coil_state(spin_system,{'Lx'},{4},'exact');

% Project out the observables
Hz=real(Hz'*traj); Hx=real(Hx'*traj);  
Nz=real(Nz'*traj); Nx=real(Nx'*traj); 
COz=real(COz'*traj); COx=real(COx'*traj);

% Plotting
kfigure(); scale_figure([1.0 2.0]);
subplot(3,1,1); plot([Hz' Hx']); 
klegend({'H$_{\rm{Z}}$','H$_{\rm{X}}$'},...
        'Location','Best'); xlim tight; kgrid;
kylabel('expectation value'); kxlabel('time, ms');
subplot(3,1,2); plot([Nz' Nx']); 
klegend({'N$_{\rm{Z}}$','N$_{\rm{X}}$'},...
        'Location','Best'); xlim tight; kgrid;
kylabel('expectation value'); kxlabel('time, ms');
subplot(3,1,3); plot([COz' COx']); 
klegend({'C$^{(\rm{O})}_{\rm{Z}}$','C$^{(\rm{O})}_{\rm{X}}$'},...
        'Location','Best'); xlim tight; kgrid;
kxlabel('time, ms'); kylabel('expectation value'); 

end

% Flip-lock-backflip-lock-backflip pumping sequence
function traj=frydman_pump(spin_system,parameters,~,R,K)

% Get pulse operators
Hy=operator(spin_system,'Ly',parameters.spins{1});
Ny=operator(spin_system,'Ly',parameters.spins{2});
Cy=operator(spin_system,'Ly',parameters.spins{3});

% Isotropic thermal equilibrium
rho=equilibrium(spin_system);

% Freeze additive replacement because all concentrations are invariant
K=K(0,unit_state(spin_system));

% Effective DIPSI Hamiltonian, CP stage 1
spin_system=dictum(spin_system,{'1H'},'ignore');         % Kill Zeeman on H
spin_system=dictum(spin_system,{'15N'},'ignore');        % Kill Zeeman on N
spin_system=dictum(spin_system,{'13C'},'ignore');        % Kill Zeeman on C
spin_system=dictum(spin_system,{'1H','15N'},'strong');   % Make H-N J-coupling strong
spin_system=dictum(spin_system,{'15N','13C'},'strong');  % Make N-C J-coupling strong
spin_system=dictum(spin_system,{'1H','13C'},'strong');   % Make H-C J-coupling strong
H_CP1=hamiltonian(spin_system); L_CP1=H_CP1+1i*R+1i*K;

% Effective DIPSI Hamiltonian, CP stage 2
spin_system=dictum(spin_system,{'1H'},'secular');        % Return Zeeman on H
spin_system=dictum(spin_system,{'1H','15N'},'ignore');   % Disconnect H-N J-coupling
spin_system=dictum(spin_system,{'1H','13C'},'ignore');   % Disconnect H-C J-coupling
H_CP2=hamiltonian(spin_system); L_CP2=H_CP2+1i*R+1i*K;

% Get the trajectory going
time_step=parameters.cp_dur./parameters.cp_npt; traj=[];

% N,C flip and crush
rho=step(spin_system,Ny,rho,+pi/2); 
rho=step(spin_system,Cy,rho,+pi/2);  
rho=homospoil(spin_system,rho,'destroy');

% Sequence loop
for n=1:parameters.nloops
    
    % Forward flips
    rho=step(spin_system,Hy,rho,+pi/2);   
    rho=step(spin_system,Cy,rho,+pi/2);   
    rho=step(spin_system,Ny,rho,+pi/2);
    
    % First CP period
    rho_stack=evolution(spin_system,L_CP1,[],rho,time_step(1),...
                        parameters.cp_npt(1)-1,'trajectory');
    traj=[traj rho_stack]; rho=rho_stack(:,end); %#ok<AGROW>
    
    % Proton back-flip
    rho=step(spin_system,Hy,rho,-pi/2);
    
    % Second CP period
    rho_stack=evolution(spin_system,L_CP2,[],rho,time_step(2),...
                        parameters.cp_npt(2)-1,'trajectory');
    traj=[traj rho_stack]; rho=rho_stack(:,end); %#ok<AGROW>
     
    % Carbon and nitrogen back-flips    
    rho=step(spin_system,Cy,rho,-pi/2);     
    rho=step(spin_system,Ny,rho,-pi/2);
    
    % Crusher gradient     
    rho=homospoil(spin_system,rho,'destroy');
    
end

end

