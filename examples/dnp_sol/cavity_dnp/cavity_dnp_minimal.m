% One ideal cavity-assisted DNP cycle of Lorriaux et al., arXiv:2610.19180v1.
% A Bi:Si donor and a Cn cavity are declared as native Spinach particles.
% The complete 20*nlevels Hilbert space is retained. The three resonant
% donor states |6>, |7>, |15> are selected spectrally; other donor levels
% remain spectators. Ideal selective pi pulses precede photon emission.
%
%                    result=cavity_dnp_minimal(nlevels)
%
% nlevels is the required cavity Fock dimension, an integer >=3. Time is
% in seconds, angular frequencies in rad/s. The result structure contains
% trajectories, native generators, observables, and consistency checks.
% This is a single-donor illustration, not an ensemble or collective fit.
% Coupling is switched off during the ideal pulses; photon loss remains.
%
% Calculation time: seconds to minutes
%
% ilya.kuprov@weizmann.ac.il

function result=cavity_dnp_minimal(nlevels)

% Check consistency
grumbler(nlevels);

% Physical donor and explicit native cavity subsystem
sys.magnet=0.06537;
sys.isotopes={'E','209Bi',['C' num2str(nlevels)]};
sys.output='hush';
sys.parallel={'processes',1};
inter.coupling.scalar=cell(3);
inter.coupling.scalar{1,2}=1.475e9;
inter.modes.frqs={[] [] 0};
inter.modes.linewidths={[] [] 210e3};
inter.temperature=0;
bas.formalism='zeeman-hilb';
bas.approximation={'none'};

% Native Spinach system and full donor Hamiltonian
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
H=hamiltonian(assume(spin_system,'labframe'));

% Diagonalise the donor vacuum block to identify the paper's energy labels
vacuum=operator(spin_system,{'BL1'},{3});
vac_idx=find(real(diag(vacuum))>0.5);
[V,D]=eig(full(H(vac_idx,vac_idx)));
[energies,order]=sort(real(diag(D))); V=V(:,order);

% Native donor observables and spectral projectors on the full composite space
Sx=operator(spin_system,'Lx',1);
Ix=operator(spin_system,'Lx',2);
Iz=operator(spin_system,'Lz',2);
a=operator(spin_system,{'A'},{3});
projectors=cell(1,20);
for n=1:20
    psi=sparse(vac_idx,ones(size(vac_idx)),V(:,n),size(H,1),1);
    projectors{n}=psi*psi';
    for k=1:nlevels-1
        psi=a'*psi/sqrt(k);
        projectors{n}=projectors{n}+psi*psi';
    end
end
active=projectors([6 7 15]);

% Select the resonant microwave and radiofrequency transitions
lower=active{1}*Sx*active{3};
mu=norm(lower,'fro')/sqrt(nlevels); lower=lower/mu;
rf_lower=active{1}*Ix*active{2};
rf_lower=rf_lower/(norm(rf_lower,'fro')/sqrt(nlevels));
Xmw=lower+lower'; Xrf=rf_lower+rf_lower';

% Native cavity operators on the full donor-cavity Hilbert space
number=operator(spin_system,{'N'},{3});
top=operator(spin_system,{['BL' num2str(nlevels)]},{3});
Hjc=lower*a'+lower'*a;

% Ground-manifold mixture and cavity vacuum built using native state
ground=sparse(size(H,1),size(H,2));
for n=1:9
    ground=ground+projectors{n};
end
rho=ground*state(spin_system,{'BL1'},{3})/9;
rho=hilb2liouv(rho,'statevec');

% Native Liouville basis, commutators, and cavity relaxation
bas.formalism='zeeman-liouv';
spin_system=basis(spin_system,bas);
R=relaxation(spin_system);
g=2*pi*100; kappa=2*pi*210e3; omega=2*pi*10e3;
Hjc=hilb2liouv(Hjc,'comm');
Lwait=g*Hjc+1i*R;
Lmw=(omega/2)*hilb2liouv(Xmw,'comm')+1i*R;
Lrf=(omega/2)*hilb2liouv(Xrf,'comm')+1i*R;

% Native time-domain propagation of the two ideal selective pulses
nsteps=32; pulse_time=pi/omega;
traj_mw=evolution(spin_system,Lmw,[],rho,pulse_time/nsteps,nsteps,'trajectory');
traj_rf=evolution(spin_system,Lrf,[],traj_mw(:,end),...
                  pulse_time/nsteps,nsteps,'trajectory');

% Native propagation with the explicit cavity retained during emission
gamma_p=4*g^2/kappa;
wait_time=linspace(0,8/gamma_p,241);
traj_wait=evolution(spin_system,Lwait,[],traj_rf(:,end),...
                    wait_time(2),numel(wait_time)-1,'trajectory');
trajectory=[traj_mw traj_rf(:,2:end) traj_wait(:,2:end)];
result.time=[linspace(0,2*pulse_time,2*nsteps+1)...
             2*pulse_time+wait_time(2:end)];

% Native observables projected against the complete composite trajectory
result.populations=zeros(3,size(trajectory,2));
for n=1:3
    result.populations(n,:)=real(active{n}(:)'*trajectory);
end
result.photons=real(number(:)'*trajectory);
result.top_population=real(top(:)'*trajectory);
result.enhancement=result.populations(1,end)/(1/9);
result.delta_iz=real(Iz(:)'*(trajectory(:,end)-trajectory(:,1)));

% Full density-matrix consistency checks
trace_error=0; hermiticity_error=0; min_eigenvalue=0;
for n=1:size(trajectory,2)
    density=reshape(trajectory(:,n),20*nlevels,20*nlevels);
    trace_error=max(trace_error,abs(trace(density)-1));
    hermiticity_error=max(hermiticity_error,norm(density-density','fro'));
    min_eigenvalue=min(min_eigenvalue,min(real(eig(full((density+density')/2)))));
end
result.checks=struct('trace_error',trace_error,...
    'hermiticity_error',hermiticity_error,'min_eigenvalue',min_eigenvalue,...
    'max_top_population',max(result.top_population));

% Exact damped Jaynes-Cummings survival for the resonant selected transition
zeta=sqrt((kappa/4)^2-g^2);
rslow=-g^2/(kappa/4+zeta); rfast=-kappa/4-zeta;
amplitude=((1+kappa/(4*zeta))*exp(rslow*wait_time)+...
           (1-kappa/(4*zeta))*exp(rfast*wait_time))/2;
result.wait_time=wait_time;
result.exact_excited=amplitude.^2/9;
result.wait_excited=result.populations(3,2*nsteps+1:end);
result.checks.exact_error=max(abs(result.wait_excited-result.exact_excited));

% Native system and generators retained for independent control calculations
result.parameters=struct('nlevels',nlevels,'g',g,'kappa',kappa);
result.transition_hz=[energies(15)-energies(6) energies(7)-energies(6)]/(2*pi);
result.electron_matrix_element=mu;
result.spin_system=spin_system; result.generator=Lwait;
result.coupling_generator=Hjc; result.relaxation=R;
result.initial_after_pulses=traj_rf(:,end);
result.final_density=reshape(traj_wait(:,end),20*nlevels,20*nlevels);
result.projectors=active; result.number=number;
result.wait_trajectory=traj_wait;

end

% Consistency enforcement
function grumbler(nlevels)
if ~isnumeric(nlevels)||~isreal(nlevels)||~isscalar(nlevels)||...
   ~isfinite(nlevels)||(nlevels<3)||(mod(nlevels,1)~=0)
    error('nlevels must be a real integer greater than or equal to three.');
end
end


