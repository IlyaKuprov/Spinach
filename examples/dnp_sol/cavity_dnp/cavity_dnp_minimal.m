% Cavity-assisted DNP of a low-field Bi:Si donor, inspired by Fig. 1
% and Appendix D of Lorriaux et al., arXiv:2610.19180v1. The full E,
% 209Bi, and Cn system remains in sphten-liouv throughout. Native
% finite RF saturation and selective pulses prepare the donor; native
% Jaynes-Cummings exchange and photon loss reset it through the cavity.
%
%           result=cavity_dnp_minimal(nlevels,max_rank,ncycles)
%
% nlevels  - cavity Fock dimension, integer at least three
% max_rank - pulse Fokker-Planck rank, default 2
% ncycles  - phase-cycled RF saturation sweeps, default 300
%
% The donor starts at 25 mK. Rectangular phase-cycled RF pulses replace
% the experimental WURST scrambling; finite rectangular selective pulses
% replace its adiabatic inversions. The cavity is disconnected during
% preparation and connected in vacuum for reset. Its bath is taken at
% zero temperature. Coupling is an illustrative single-donor value, not
% an ensemble fit. Nuclear drive through its much smaller magnetic moment
% is neglected. The complete isotropic hyperfine interaction is retained.
%
% Calculation time: minutes; requests a 64-worker parallel pool
%
% ilya.kuprov@weizmann.ac.il

function result=cavity_dnp_minimal(nlevels,max_rank,ncycles)

% Defaults and consistency checks
if nargin<2, max_rank=2; end
if nargin<3, ncycles=300; end
grumbler(nlevels,max_rank,ncycles);

% Actual low-field donor and an initially disconnected bosonic spectator
sys.magnet=0.06537; sys.isotopes={'E','209Bi',['C' num2str(nlevels)]};
sys.output='hush'; sys.parallel={'local',64}; sys.enable={'paranoia'};
inter.zeeman.scalar={2.00032,0,[]};
inter.coupling.scalar=cell(3,3); inter.coupling.scalar{1,2}=1.475e9;
inter.modes.frqs={[] [] 0}; inter.temperature=0.025;
bas.formalism='sphten-liouv'; bas.approximation={'none'};

% Native laboratory-frame Hamiltonian and thermal donor preparation
s=assume(basis(create(sys,inter),bas),'labframe');
H=hamiltonian(s); Hl=hamiltonian(s,'left'); rho=equilibrium(s,Hl);
unit=state(s,'E',1); Lx=operator(s,'Lx',1); Ly=operator(s,'Ly',1);

% Exact native invariant-space reduction removes the disconnected mode
blocks=reduce(s,H+Lx+Ly,rho); P=horzcat(blocks{:});
H=P'*H*P; Hl=P'*Hl*P; rho=P'*rho;
Lx=P'*Lx*P; Ly=P'*Ly*P; unit_reduced=P'*unit;
Iz=P'*operator(s,'Lz',2,'left')'*unit;

% Native left-product eigenvalues supply frequencies, not eigenstates
energies=sort(real(eig(full(Hl))));
levels=uniquetol(energies,size(Hl,1)*eps(norm(Hl,1)),DataScale=1);
assert(numel(levels)==20,'The full donor must have twenty distinct levels.');
rf_frqs=diff(levels(1:9))/(2*pi);
mw_frq=(levels(15)-levels(6))/(2*pi);
rf_frq=(levels(7)-levels(6))/(2*pi);

% Native finite-pulse maps for RF saturation of all eight ground lines
rf_maps=cell(8,1); rf_amplitude=2*pi*100e3;
parfor line=1:8
    [~,~,rf_maps{line}]=shaped_pulse_af(s,H,Lx,0*Ly,rho,...
        rf_frqs(line),rf_amplitude,10e-6,0,max_rank,'expm');
end

% Native coherence selection implements a uniform donor RF phase ensemble
for cycle=1:ncycles
    previous=rho;
    for line=1:8
        rho=P'*coherence(s,P*(rf_maps{line}*rho),{{[1 2],0}});
    end
end
scramble_change=norm(rho-previous);
preparation_trace=real(unit_reduced'*rho);
rho=rho/(unit_reduced'*rho); initial_iz=real(Iz'*rho);

% Finite selective microwave pulse, including the full donor drift
mw_time=17.02542402393068e-6;
rho_mw=shaped_pulse_af(s,H,Lx,0*Ly,rho,mw_frq,...
    rf_amplitude,mw_time,0,max_rank,'expm');
rho_mw=P'*coherence(s,P*rho_mw,{{[1 2],0}});

% Finite selective RF pulse, with the same native phase-ensemble selection
rf_time=22.56805013451008e-6;
rho_rf=shaped_pulse_af(s,H,Lx,0*Ly,rho_mw,rf_frq,...
    rf_amplitude,rf_time,0,max_rank,'expm');
rho_rf=P'*coherence(s,P*rho_rf,{{[1 2],0}});
pulse_iz=real(Iz'*[rho_mw rho_rf]);

% Attach vacuum with native bosonic syntax; retain controls for validation
vacuum=operator(s,'BL1',3,'left');
rho_wait=vacuum*(P*rho_rf); rho_no_rf=vacuum*(P*rho_mw);
rho_initial=vacuum*(P*rho);
rho_wait=rho_wait/(unit'*rho_wait);
rho_no_rf=rho_no_rf/(unit'*rho_no_rf);
rho_initial=rho_initial/(unit'*rho_initial);

% Rebuild the connected system from native interaction specifications
inter.modes.frqs={[] [] mw_frq};
inter.modes.exchange=cell(3,3); inter.modes.exchange{1,3}=100;
inter.modes.linewidths={[] [] 210e3}; inter.temperature=0;
s=assume(basis(create(sys,inter),bas),'labframe');
s.inter.modes.strength.exchange='rwa';
H=hamiltonian(s); R=relaxation(s); unit=state(s,'E',1);

% Select the native MATLAB-expm branch for these bounded reset matrices
s.tols.small_matrix=size(H,1)+1;

% Native physical trace duals, with no Euclidean-metric assumption
Iz=operator(s,'Lz',2,'left')'*unit;
number=operator(s,'N',3,'left')'*unit;
top=operator(s,['BL' num2str(nlevels)],3,'left')'*unit;
observables=[unit Iz number top];

% Laboratory-frame reset: retain GHz donor energies and explicit photons
time=0:0.001:20;
signals=evolution(s,H+1i*R,observables,rho_wait,time(2),...
                  numel(time)-1,'multichannel');
signals=real(signals);

% Return the measured trajectories and native objects for independent tests
result=struct('time',time,'trace',signals(1,:),'iz',signals(2,:),...
    'photons',signals(3,:),'top_population',signals(4,:),...
    'initial_iz',full(initial_iz),'pulse_iz',full(pulse_iz),...
    'delta_iz',full(signals(2,end)-initial_iz),...
    'scramble_change',scramble_change,'preparation_trace',full(preparation_trace),...
    'transition_hz',[mw_frq rf_frq],'saturation_hz',rf_frqs,...
    'pulse_durations',[mw_time rf_time],'spin_system',s,...
    'H',H,'R',R,'rho',rho_wait,'rho_no_rf',rho_no_rf,...
    'rho_initial',rho_initial,'observables',observables,...
    'nlevels',nlevels,'max_rank',max_rank,'ncycles',ncycles);

end

% Consistency enforcement
function grumbler(nlevels,max_rank,ncycles)
if ~isnumeric(nlevels)||~isreal(nlevels)||~isscalar(nlevels)||...
   ~isfinite(nlevels)||(nlevels<3)||(mod(nlevels,1)~=0)
    error('nlevels must be an integer greater than or equal to three.');
end
if ~isnumeric(max_rank)||~isreal(max_rank)||~isscalar(max_rank)||...
   ~isfinite(max_rank)||(max_rank<1)||(mod(max_rank,1)~=0)
    error('max_rank must be a positive integer.');
end
if ~isnumeric(ncycles)||~isreal(ncycles)||~isscalar(ncycles)||...
   ~isfinite(ncycles)||(ncycles<1)||(mod(ncycles,1)~=0)
    error('ncycles must be a positive integer.');
end
end

