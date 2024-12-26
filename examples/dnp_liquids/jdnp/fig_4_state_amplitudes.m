% Time evolution of the individual states in the basis set built
% from the singlet-triplet basis on the two electrons and Carte-
% sian spin operator basis on the nucleus, demonstrating the im-
% balance between singlet-alpha and singlet-beta subspace. This
% leads to trainsient nuclear spin polarisation. Details in:
%
%               https://doi.org/10.1039/d1cp04186j
%
% Calculation time: seconds
%
% mariagrazia.concilio@sjtu.edu.cn
% ilya.kuprov@weizmann.ac.il

function fig_4_state_amplitudes()

% Load the spin system
[sys,inter,bas,parameters]=system_specification();

% Experiment parameters
parameters.mw_pwr=2*pi*250e3;
parameters.t_step=1e-3;
parameters.nsteps=700;

% Set the magnet
sys.magnet=14.08;

% Set microwave offset frequency
f_free=g2freq(parameters.g_ref,sys.magnet);
f_trityl=g2freq(parameters.g_trityl,sys.magnet);
parameters.mw_off=2*pi*(f_trityl-f_free);   
        
% Set the exchange coupling
electron_zeeman_tensor=g2freq(inter.zeeman.matrix{2},sys.magnet);
electron_zeeman_iso=mean(diag(electron_zeeman_tensor));
proton_zeeman_iso=sys.magnet*spin('1H')/(2*pi);
inter.coupling.scalar{2,3}=electron_zeeman_iso+proton_zeeman_iso;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Electron operators
Ex=operator(spin_system,'Lx','E');
Ez=operator(spin_system,'Lz','E');

% Thermal equilibirium state
H0=hamiltonian(assume(spin_system,'labframe'),'left');
rho_eq=equilibrium(spin_system,H0);

% Hamiltonian and relaxation superoperator
H=hamiltonian(assume(spin_system,'esr'));
R=relaxation(spin_system);

% Add microwave terms to the Hamiltonian
H=H+parameters.mw_pwr*Ex+parameters.mw_off*Ez;
    
% Build component operators
unit=unit_state(spin_system);
Nz=state(spin_system,{'Lz'},{1});
E1z=state(spin_system,{'Lz'},{2});
E2z=state(spin_system,{'Lz'},{3});
NzE1z=state(spin_system,{'Lz','Lz'},{1,2});
NzE2z=state(spin_system,{'Lz','Lz'},{1,3});
LzE2z=state(spin_system,{'Lz','Lz'},{2,3});
E1mE2p=state(spin_system,{'L-','L+'},{2,3});
E1pE2m=state(spin_system,{'L+','L-'},{2,3});
NzE1zE2z=state(spin_system,{'Lz','Lz','Lz'},{1,2,3});
NzE1mE2p=state(spin_system,{'Lz','L-','L+'},{1,2,3});
NzE1pE2m=state(spin_system,{'Lz','L+','L-'},{1,2,3});

% Build alpha singlet and triplet states
Sa=(unit/8)+(Nz/4)-(E1mE2p/4)-(E1pE2m/4)-(LzE2z/2)-(NzE1mE2p/2)-(NzE1pE2m/2)-NzE1zE2z;
Tpa=(unit/8)+(Nz/4)+(E1z/4)+(E2z/4)+(LzE2z/2)+(NzE2z/2)+(NzE1z/2)+NzE1zE2z;
T0a=(unit/8)+(Nz/4)+(E1mE2p/4)+(E1pE2m/4)-(LzE2z/2)+(NzE1mE2p/2)+(NzE1pE2m/2)-NzE1zE2z;
Tma=(unit/8)+(Nz/4)-(E1z/4)-(E2z/4)+(LzE2z/2)-(NzE2z/2)-(NzE1z/2)+NzE1zE2z;

% Build beta singlet and triplet states
Sb=(unit/8)-(Nz/4)-(E1mE2p/4)-(E1pE2m/4)-(LzE2z/2)+(NzE1mE2p/2)+(NzE1pE2m/2)+NzE1zE2z;
Tpb=(unit/8)-(Nz/4)+(E1z/4)+(E2z/4)+(LzE2z/2)-(NzE2z/2)-(NzE1z/2)-NzE1zE2z;
T0b=(unit/8)-(Nz/4)+(E1mE2p/4)+(E1pE2m/4)-(LzE2z/2)-(NzE1mE2p/2)-(NzE1pE2m/2)+NzE1zE2z;
Tmb=(unit/8)-(Nz/4)-(E1z/4)-(E2z/4)+(LzE2z/2)+(NzE2z/2)+(NzE1z/2)-NzE1zE2z;
    
% Build Nz-singlet and Nz-triplet states
SNz=(Nz/4)-(NzE1mE2p/2)-(NzE1pE2m/2)-NzE1zE2z;
TpNz=(Nz/4)+(NzE1z/2)+(NzE2z/2)+NzE1zE2z;
T0Nz=(Nz/4)+(NzE1mE2p/2)+(NzE1pE2m/2)-NzE1zE2z;
TmNz=(Nz/4)-(NzE1z/2)-(NzE2z/2)+NzE1zE2z;

% Detection states
coils=[Tpa,Tpb,Tma,Tmb,T0a,T0b,Sa,Sb,SNz,TpNz,T0Nz,TmNz,E1z,E2z,Nz];
    
% Run the time evolution
answer=evolution(spin_system,H+1i*R,coils,rho_eq,parameters.t_step,...
                 parameters.nsteps,'multichannel');

% Time axis generation
t_axis=linspace(0,parameters.t_step*parameters.nsteps,parameters.nsteps+1);   
                 
% Plotting
figure(); subplot(1,3,1); scale_figure([1.75 0.5]);
plot(t_axis(2:end),real(answer(1,2:end)),'b-'); hold on;
plot(t_axis(2:end),real(answer(3,2:end)),'b-');
plot(t_axis(2:end),real(answer(5,2:end)),'b-');
plot(t_axis(2:end),real(answer(2,2:end)),'r-');
plot(t_axis(2:end),real(answer(4,2:end)),'r-');
plot(t_axis(2:end),real(answer(6,2:end)),'r-'); kgrid;
xlabel('time / seconds','interpreter','latex');
ylabel('state population','interpreter','latex');
legend({'${T_{+,\alpha}}$','${T_{-,\alpha}}$','${T_{0,\alpha}}$',...
        '${T_{+,\beta}}$','${T_{-,\beta}}$','${T_{0,\beta}}$'},...
        'interpreter','latex','Location','northeast');
set(gca,'TickLabelInterpreter','latex'); axis tight;
subplot(1,3,2);
plot(t_axis(2:end),real(answer(7,2:end)),'b-'); hold on;
plot(t_axis(2:end),real(answer(8,2:end)),'r-'); kgrid;
xlabel('time / seconds','interpreter','latex');
ylabel('state population','interpreter','latex');
legend({'${S_{\alpha}}$','${S_{\beta}}$'},...
        'interpreter','latex','Location','southeast');
set(gca,'TickLabelInterpreter','latex'); axis tight;
subplot(1,3,3);
plot(t_axis(2:end),real(answer(15,2:end)),'b-'); hold on;
xlabel('time / seconds','interpreter','latex');
ylabel('state population','interpreter','latex');
legend({'${N_{\rm{Z}}}$'},'interpreter','latex','Location','southeast');
set(gca,'TickLabelInterpreter','latex'); kgrid; axis tight;

end

