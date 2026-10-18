% Minimal cavity-assisted DNP cycle of Lorriaux et al., arXiv:2610.19180v1.
% Three active Bi:Si eigenstates |6>, |7>, |15>, plus seven inert ground
% populations. The cavity is a finite Weyl boson, not a Purcell-rate term.
% Ideal selective pi pulses replace the experimental adiabatic pulses.
% Call after adding Spinach to the MATLAB path:
%
%                     result=cavity_dnp_minimal(nlevels)
%
% nlevels is the required number of cavity Fock states (integer >=2).
% This single-donor illustration omits ensemble and collective effects.
% The cavity coupling is switched off during the ideal selective pulses.
% Rates are rad/s or inverse seconds; time is seconds. The result structure
% contains the trajectory and diagnostics. No figures or files are created.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function result=cavity_dnp_minimal(nlevels)

% Check consistency
grumbler(nlevels);

% Full Bi:Si Hamiltonian for the physical state identification
sys.magnet=0;
sys.isotopes={'E','209Bi'};
sys.output='hush';
inter=struct();
bas.formalism='zeeman-hilb';
bas.approximation={'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Electron and nuclear Cartesian operators
Sx=operator(spin_system,'Lx',1); Ix=operator(spin_system,'Lx',2);
Sy=operator(spin_system,'Ly',1); Iy=operator(spin_system,'Ly',2);
Sz=operator(spin_system,'Lz',1); Iz=operator(spin_system,'Lz',2);
Hbi=2*pi*(27.997e9*0.06537*Sz-6.962e6*0.06537*Iz+...
          1.475e9*(Sx*Ix+Sy*Iy+Sz*Iz));
[V,E]=eig(full(Hbi),'vector');
[energies,order]=sort(real(E)); V=V(:,order);

% Active eigenstates in the order |6>, |7>, |15>
active=V(:,[6 7 15]);
result.transition_hz=[energies(15)-energies(6);...
                      energies(7)-energies(6)]/(2*pi);
result.active_iz=real(diag(active'*Iz*active));
result.electron_matrix_element=abs(V(:,6)'*Sx*V(:,15));

% Explicit boson and selective transition matrices
cavity=weyl(nlevels);
a=kron(speye(3),cavity.a);
number=a'*a;
lower=sparse(1,3,1,3,3);
Hjc=kron(lower,cavity.c)+kron(lower',cavity.a);
Xmw=kron(lower+lower',cavity.u);
Xrf=kron(sparse([1 2],[2 1],[1 1],3,3),cavity.u);

% Demonstration parameters: measured cavity width, illustrative one-spin g
kappa=2*pi*210e3;
g=2*pi*100;
omega_mw=2*pi*10e3;
omega_rf=2*pi*10e3;
gamma_p=4*g^2/kappa;
result.parameters=struct('nlevels',nlevels,'kappa',kappa,'g',g,...
    'omega_mw',omega_mw,'omega_rf',omega_rf,'gamma_p',gamma_p);

% Vectorised Lindblad generator: column-major vec convention
identity=speye(3*nlevels);
comm=@(H) kron(identity,H)-kron(H.',identity);
R=kappa*(kron(conj(a),a)-0.5*kron(identity,number)-...
         0.5*kron(number.',identity));
Lwait=g*comm(Hjc)+1i*R;
Lmw=0.5*omega_mw*comm(Xmw)+1i*R;
Lrf=0.5*omega_rf*comm(Xrf)+1i*R;

% Initially equal ground populations: active block trace 2/9, spectators 7/9
vacuum=sparse(1,1,1,nlevels,nlevels);
rho=kron(diag([1 1 0])/9,vacuum);
rho=rho(:);

% Time-resolved microwave and radio-frequency pi pulses
nsteps=32;
tmw=pi/omega_mw; trf=pi/omega_rf;
Pmw=propagator(spin_system,Lmw,tmw/nsteps);
Prf=propagator(spin_system,Lrf,trf/nsteps);
trajectory=zeros(numel(rho),2*nsteps+1); trajectory(:,1)=rho;
for n=1:nsteps
    rho=Pmw*rho; trajectory(:,n+1)=rho;
end
for n=1:nsteps
    rho=Prf*rho; trajectory(:,nsteps+n+1)=rho;
end
pulse_time=[linspace(0,tmw,nsteps+1) tmw+(1:nsteps)*trf/nsteps];

% Wait for eight nominal Purcell lifetimes; cavity remains in the model
wait_time=linspace(0,8/gamma_p,241);
Pwait=propagator(spin_system,Lwait,wait_time(2));
wait_trajectory=zeros(numel(rho),numel(wait_time));
wait_trajectory(:,1)=rho;
for n=2:numel(wait_time)
    rho=Pwait*rho; wait_trajectory(:,n)=rho;
end
result.time=[pulse_time tmw+trf+wait_time(2:end)];
trajectory=[trajectory wait_trajectory(:,2:end)];

% Populations, photons, and physical density-matrix checks
projectors=cell(1,3);
for n=1:3
    projectors{n}=kron(sparse(n,n,1,3,3),cavity.u);
end
top=kron(speye(3),sparse(nlevels,nlevels,1,nlevels,nlevels));
result.populations=zeros(3,size(trajectory,2));
result.photons=zeros(1,size(trajectory,2));
result.top_population=zeros(1,size(trajectory,2));
trace_error=0; hermiticity_error=0; min_eigenvalue=0;
for n=1:size(trajectory,2)
    density=reshape(trajectory(:,n),3*nlevels,3*nlevels);
    for k=1:3
        result.populations(k,n)=real(trace(projectors{k}*density));
    end
    result.photons(n)=real(trace(number*density));
    result.top_population(n)=real(trace(top*density));
    trace_error=max(trace_error,abs(trace(density)-2/9));
    hermiticity_error=max(hermiticity_error,norm(density-density','fro'));
    min_eigenvalue=min(min_eigenvalue,min(real(eig((density+density')/2))));
end
result.enhancement=9*result.populations(1,end);
result.delta_iz=result.active_iz.'*(result.populations(:,end)-...
                                   result.populations(:,1));
result.checks=struct('trace_error',trace_error,...
    'hermiticity_error',hermiticity_error,'min_eigenvalue',min_eigenvalue,...
    'max_top_population',max(result.top_population));

% Closed-form damped Jaynes-Cummings survival, no cavity elimination
zeta=sqrt((kappa/4)^2-g^2);
rslow=-g^2/(kappa/4+zeta); rfast=-kappa/4-zeta;
amplitude=((1+kappa/(4*zeta))*exp(rslow*wait_time)+...
           (1-kappa/(4*zeta))*exp(rfast*wait_time))/2;
result.wait_time=wait_time;
result.exact_excited=amplitude.^2/9;
result.purcell_excited=exp(-gamma_p*wait_time)/9;
result.wait_excited=result.populations(3,2*nsteps+1:end);
result.checks.exact_error=max(abs(result.wait_excited-result.exact_excited));
result.checks.purcell_error=max(abs(result.wait_excited-result.purcell_excited));
result.final_density=reshape(rho,3*nlevels,3*nlevels);
result.generator=Lwait;
result.initial_after_pulses=wait_trajectory(:,1);
result.spin_system=spin_system;

end

% Consistency enforcement
function grumbler(nlevels)
if ~isnumeric(nlevels)||~isreal(nlevels)||~isscalar(nlevels)||...
   ~isfinite(nlevels)||(nlevels<2)||(mod(nlevels,1)~=0)
    error('nlevels must be a real integer greater than or equal to two.');
end
end


