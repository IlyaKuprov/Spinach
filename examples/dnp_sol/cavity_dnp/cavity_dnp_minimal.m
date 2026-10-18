% Minimal cavity-assisted nuclear pumping cycle inspired by Lorriaux
% et al., arXiv:2610.19180v1, Fig. 1. An electron, a spin-1/2 nucleus,
% and an explicit damped bosonic cavity illustrate population shelving,
% a nuclear swap, and radiative reset. This is a high-field four-level
% analogue, NOT the low-field 20-level Bi:Si donor or an ensemble fit.
%
%                    result=cavity_dnp_minimal(nlevels)
%
% nlevels is the cavity Fock dimension, an integer at least three. All
% states, generators, pulses, relaxation, and evolution use the native
% sphten-liouv formalism. Frequencies are in Hz, times in seconds.
% Ideal conditional pulses omit drift and decay during the pulses.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function result=cavity_dnp_minimal(nlevels)

% Check the cavity truncation
grumbler(nlevels);

% An effective electron-nuclear doublet and a quantised cavity
sys.magnet=0; sys.isotopes={'E','1H',['C' num2str(nlevels)]};
sys.output='hush'; sys.parallel={'local',2}; sys.enable={'paranoia'};
bas.formalism='sphten-liouv'; bas.approximation={'none'};

% Common rotating frame: secular hyperfine splitting and resonant mode
inter.coupling.scalar=cell(3,3); inter.coupling.scalar{1,2}=10e3;
inter.modes.frqs={[] [] 5e3}; inter.modes.linewidths={[] [] 2e3};
inter.modes.exchange=cell(3,3); inter.modes.exchange{1,3}=100;
inter.temperature=0;

% Native interaction specification generates the entire Liouvillian
spin_system=create(sys,inter); spin_system=basis(spin_system,bas);
spin_system=assume(spin_system,'cavity');
spin_system.inter.coupling.strength{1,2}='zz';
H=hamiltonian(spin_system); R=relaxation(spin_system);

% Ground electron with unpolarised nucleus and vacuum cavity
rho=state(spin_system,{'ZL1','BL1'},{1,3});
unit=state(spin_system,'E','all'); rho=rho/(unit'*rho);

% Native observables: ground electron with nuclear alpha and beta
target=operator(spin_system,{'ZL1','ZL2'},{1,2},'left')'*unit;
source=operator(spin_system,{'ZL1','ZL1'},{1,2},'left')'*unit;
excited=operator(spin_system,{'ZL2','ZL2'},{1,2},'left')'*unit;
iz=operator(spin_system,'Lz',2,'left')'*unit;
number=operator(spin_system,'N',3,'left')'*unit;
top=operator(spin_system,['BL' num2str(nlevels)],3,'left')'*unit;
observables=[target source excited iz number top unit];
initial=real(observables'*rho);

% Ideal microwave pi pulse on the nuclear-alpha electron transition
microwave=operator(spin_system,{'Lx','ZL2'},{1,2});
rho=step(spin_system,microwave,rho,pi);

% Ideal radiofrequency pi pulse within the electron-ground doublet
radiofrequency=operator(spin_system,{'ZL1','Lx'},{1,2});
rho=step(spin_system,radiofrequency,rho,pi);
after_pulses=real(observables'*rho);

% Explicit cavity exchange and photon loss during the reset interval
time=linspace(0,0.1,401);
trajectory=evolution(spin_system,H+1i*R,[],rho,time(2),...
                     numel(time)-1,'trajectory');
signals=real(observables'*trajectory);

% Return observables and native objects for independent control runs
result=struct('time',time,'populations',full(signals(1:3,:)),...
              'iz',full(signals(4,:)),'photons',full(signals(5,:)),...
              'top_population',full(signals(6,:)),'trace',full(signals(7,:)),...
              'initial',full(initial),'after_pulses',full(after_pulses),...
              'enhancement',full(signals(1,end)/initial(1)),...
              'delta_iz',full(signals(4,end)-initial(4)),...
              'spin_system',spin_system,'H',H,'R',R,'rho',rho,...
              'observables',observables,'nlevels',nlevels);

end

% Consistency enforcement
function grumbler(nlevels)
if ~isnumeric(nlevels)||~isreal(nlevels)||~isscalar(nlevels)||...
   ~isfinite(nlevels)||(nlevels<3)||(mod(nlevels,1)~=0)
    error('nlevels must be a real integer greater than or equal to three.');
end
end

