% Ideal Hahn echo; one pi pulse at half the elapsed time.
% Syntax:
%
%           [E,aux]=hahn_echo(spin_system,parameters,H,R,K)
%
% parameters.times contains total sequence durations in seconds. transition
% selects the observed two-state transition; rho0/coil or preparation may
% override the prepared density and detector. normalized defaults to true
% and divides by the actual detected initial reference. nloops is required
% for CPMG. bath_pumped_hahn requires pump.center, halfwidth, units (Hz or
% rad/s), fraction and allow_empty; optional sites/species restrict the bath.
% Ideal refocusing acts on the selected transition about x. H,R,K are the
% conventional generators or explicitly tagged projected branch generators.
%
% output=density returns column-vectorised generalized observed density.
% It requires a physical initial density and does not normalize a receiver.
%
% Optional aux contains spectral/error diagnostics for coherent projected echoes.
%
% ilya.kuprov@weizmann.ac.il

function [E,aux]=hahn_echo(spin_system,parameters,H,R,K)

% Preparation and detection are owned by this experiment; projected states are coherence blocks
grumble(spin_system,parameters,H,R,K);
if ~isfield(parameters,'normalized'), parameters.normalized=true; end
density=isfield(parameters,'output')&&strcmp(parameters.output,'density');
if density
    parameters.normalized=false;
    if ~isfield(parameters,'coil'), parameters.coil=eye(spin_system.observed_dimension); end
end
[rho0,coil,reference,meta]=cce_observation(spin_system,parameters);
if ~iscell(coil), coil={coil}; end
projected=isstruct(H); E=complex(zeros(numel(parameters.times),numel(reference)));
aux=struct();
if projected&&~density&&isempty(H.dissipator)&&~isfield(parameters,'pump')&&...
   all(isfield(H,{'hamiltonian_abs_error_bound','hamiltonian_norm'}))
    [E,aux]=projected_echo(spin_system,parameters,H,reference,meta); return
end
if density
    if projected, error('hahn_echo:density','Projected blocks do not define an observed density.'); end
    E=complex(zeros(numel(parameters.times),spin_system.observed_dimension^2));
end
if projected, generator=H; else, generator=H+1i*R+1i*K; end

% Build ideal observed-spin transition pulses without a global bath basis
model=spin_system.model; dc=spin_system.observed_dimension;
db=spin_system.bath_dimension;
if ~isfield(model,'transition')
    error('hahn_echo:transition','An explicit transition is required for ideal refocusing.');
end
a=model.transition.alpha; b=model.transition.beta;
rotation=eye(dc)-a*a'-b*b'-1i*(a*b'+b*a');
rotation=kron(rotation,eye(db));
fractions=0.5;

% At zero elapsed time no control is applied, matching the reference signal
for n=1:numel(parameters.times)
    rho=rho0; previous=0; swapped=false;
    if parameters.times(n)>0
        for j=1:numel(fractions)
            segment=generator;
            if projected&&swapped
                segment.left=H.right; segment.right=H.left;
                if ~isempty(H.dissipator)
                    segment.dissipator=H.swap_matrix*conj(H.dissipator)*H.swap_matrix';
                end
            end
            duration=(fractions(j)-previous)*parameters.times(n);
            rho=cce_evolution(spin_system,segment,[],rho,duration,1,'final');
            if projected
                swapped=~swapped;
            else
                rho=reshape(rho,dc*db,dc*db); rho=rotation*rho*rotation';
                if strcmp(spin_system.formalism,'zeeman-liouv'), rho=rho(:); end
            end
            previous=fractions(j);
        end
    end
    segment=generator;
    if projected&&swapped
        segment.left=H.right; segment.right=H.left;
        if ~isempty(H.dissipator)
            segment.dissipator=H.swap_matrix*conj(H.dissipator)*H.swap_matrix';
        end
    end
    rho=cce_evolution(spin_system,segment,[],rho,(1-previous)*parameters.times(n),1,'final');
    if density
        rho=reshape(rho,dc*db,dc*db); observed=complex(zeros(dc));
        for a=1:dc
            for b=1:dc
                observed(a,b)=trace(rho((a-1)*db+(1:db),(b-1)*db+(1:db)));
            end
        end
        E(n,:)=reshape(observed,1,[]); continue
    end
    for k=1:numel(coil)
        if isvector(rho), value=coil{k}'*rho;
        else, value=trace(coil{k}'*rho); end
        if meta.normalized, value=value/reference(k); end
        E(n,k)=value;
    end
end

end

% Validate the registered sequence contract before preparation
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(spin_system)||~isfield(spin_system,'kind')||...
   ~strcmp(spin_system.kind,'cce_subsystem')||~isstruct(parameters)||...
   ~isfield(parameters,'times')||~isnumeric(parameters.times)||...
   ~isreal(parameters.times)||~isvector(parameters.times)||...
   any(~isfinite(parameters.times))||any(parameters.times<0)||any(diff(parameters.times)<0)
    error('hahn_echo:input','A native retained subsystem and increasing nonnegative times are required.');
end
if isfield(parameters,'sequence')&&~isempty(parameters.sequence)
    error('hahn_echo:sequence','This named experiment defines its own controls.');
end
if ~(isnumeric(H)||isstruct(H))||~isnumeric(R)||~isnumeric(K)
    error('hahn_echo:generators','Supply conventional or tagged projected generators.');
end

end

% Propagate the two invariant branches once per interval; retain spectral bounds
function [E,aux]=projected_echo(system,parameters,H,reference,meta)

% Resolve the branch spectra once, not once per observation time
[va,ea]=eig(full((H.left+H.left')/2),'vector');
[vb,eb]=eig(full((H.right+H.right')/2),'vector');
ra=norm(H.left*va-va.*ea.','fro'); rb=norm(H.right*vb-vb.*eb.','fro');
oa=norm(va'*va-eye(size(H.left)),'fro'); ob=norm(vb'*vb-eye(size(H.right)),'fro');
residual=max(ra,rb); orthogonality=max(oa,ob); db=system.bath_dimension;

% Build the physical ideal pulse in the selected observed transition
transition=system.model.transition;
generator=(transition.alpha*transition.beta'+transition.beta*transition.alpha')/2;
pulse=expm(-1i*(pi*kron(full(generator),eye(db))));
pa=kron(transition.alpha,eye(db)); pb=kron(transition.beta,eye(db));
pair={pa'*pulse*pb,pb'*pulse*pa};
rho_b=1;
for id=system.sites
    where=find(system.configuration.bath_ids==id,1);
    rho_b=kron(rho_b,system.configuration.bath_rho{where});
end

% Keep the branch-product order and zero-time convention explicit
E=complex(zeros(numel(parameters.times),1)); bounds=zeros(size(E));
for n=1:numel(parameters.times)
    time=parameters.times(n);
    if time==0
        value=1;
    else
        ua=eye(db); ub=eye(db); dt=time/2;
        ua=((va.*exp(-1i*ea.'*dt))*va')*ua;
        ub=((vb.*exp(-1i*eb.'*dt))*vb')*ub;
        ua=pair{1}*ua; ub=pair{2}*ub;
        ua=((vb.*exp(-1i*eb.'*dt))*vb')*ua;
        ub=((va.*exp(-1i*ea.'*dt))*va')*ub;
        value=trace(ua*rho_b*ub');
    end
    E(n)=value;
    operations=20*db^3; gamma=operations*eps/(1-operations*eps);
    bounds(n)=2*time*(residual+H.hamiltonian_abs_error_bound)+8*orthogonality+...
              gamma*(1+time*H.hamiltonian_norm);
end

% A selected-transition detector scales only unnormalised responses
if ~meta.normalized
    E=E*reference; bounds=bounds*abs(reference);
end
aux=struct('abs_error_bound',bounds,'eigen_residual',residual,...
           'orthogonality_residual',orthogonality,...
           'hamiltonian_bound',H.hamiltonian_abs_error_bound,...
           'backend',parameters.cce_aux.backend);

end


