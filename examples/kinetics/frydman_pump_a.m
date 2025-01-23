% Lucio Frydman's water exchange based spin-lock pump, Figure 2
% from https://doi.org/10.1016/j.jmr.2021.107083
%
% Calculation time: seconds.
%
% mihajlo.novakovic@weizmann.ac.il
% lucio.frydman@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il
% mariagrazia.concilio@sjtu.edu.cn

function frydman_pump_a()

% Number of water protons
n_water_protons=20;

% Magnet field
sys.magnet=11.7; 

% Core spin system 
sys.isotopes={'1H','15N','13C','13C'};

% Add water protons
sys.isotopes(5:(5+n_water_protons-1))={'1H'};

% Chemical shifts
inter.zeeman.scalar={0.0,0.0,0.0,0.0};
inter.zeeman.scalar(5:(5+n_water_protons-1))={0.0};

% Scalar couplings, Hz
inter.coupling.scalar=cell(n_water_protons+4);
inter.coupling.scalar{1,2}=-45;  

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

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.level=4;
bas.space_level=1;

% Exchange rates, Hz
nh_wt_exch_rate=10;  % between NH and nearest water
wt_wt_exch_rate=1e4; % between all water protons

% Exchange rate matrix
inter.chem.flux_rate=zeros(n_water_protons+4);
inter.chem.flux_rate(1,5)=nh_wt_exch_rate;
inter.chem.flux_rate(5,1)=nh_wt_exch_rate;
inter.chem.flux_rate(5:(5+n_water_protons-1),...
                     5:(5+n_water_protons-1))=wt_wt_exch_rate;
inter.chem.flux_type='intermolecular';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation parameters
parameters.spins={'1H','15N','13C'};
parameters.cp_dur=100e-3;
parameters.cp_npt=100;

% Get system trajectory (sequence is below)
traj=liquid(spin_system,@frydman_pump,parameters,'nmr');

% Observables: peptide bond N-H
Hz=state(spin_system,{'Lz'},{1});
Hx=state(spin_system,{'Lx'},{1}); 
Nz=state(spin_system,{'Lz'},{2});
Nx=state(spin_system,{'Lx'},{2});

% Project out the observables
Hz=real(Hz'*traj); Hx=real(Hx'*traj);  
Nz=real(Nz'*traj); Nx=real(Nx'*traj); 

% Plotting
figure(); scale_figure([1.0 1.5]);
subplot(2,1,1); plot([Hz' Hx']); 
klegend({'H$_{\rm{Z}}$','H$_{\rm{X}}$'}); xlim tight; kgrid;
kylabel('expectation value'); kxlabel('time, ms');
subplot(2,1,2); plot([Nz' Nx']); 
klegend({'N$_{\rm{Z}}$','N$_{\rm{X}}$'}); xlim tight; kgrid;
kylabel('expectation value'); kxlabel('time, ms');

end

% Flip-and-lock basic version of the pump
function traj=frydman_pump(spin_system,parameters,~,R,K)

% Get pulse operators
Hy=operator(spin_system,'Ly',parameters.spins{1});
Ny=operator(spin_system,'Ly',parameters.spins{2});

% Lab frame Hamiltonian and equilibrium state
H=hamiltonian(assume(spin_system,'labframe'),'left');
rho=equilibrium(spin_system,H);

% Effective spin-lock Hamiltonian
spin_system=dictum(spin_system,{'1H'},'ignore');         % Kill Zeeman on H
spin_system=dictum(spin_system,{'15N'},'ignore');        % Kill Zeeman on N
spin_system=dictum(spin_system,{'13C'},'ignore');        % Kill Zeeman on C
spin_system=dictum(spin_system,{'1H','15N'},'strong');   % Make H-N J-coupling strong
spin_system=dictum(spin_system,{'15N','13C'},'strong');  % Make N-C J-coupling strong
spin_system=dictum(spin_system,{'1H','13C'},'strong');   % Make H-C J-coupling strong
H_CP1=hamiltonian(spin_system); L_CP1=H_CP1+1i*R+1i*K;

% Get the trajectory going
time_step=parameters.cp_dur/parameters.cp_npt;

% Nitrogen crusher
rho=step(spin_system,Ny,rho,+pi/2);  
rho=homospoil(spin_system,rho,'destroy');
   
% Forward flip
rho=step(spin_system,Hy,rho,+pi/2);
rho=step(spin_system,Ny,rho,+pi/2);

% CP period
traj=krylov(spin_system,L_CP1,[],rho,time_step,...
            parameters.cp_npt-1,'trajectory');

end

