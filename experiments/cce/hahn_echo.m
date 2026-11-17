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
% Generalized unitary Zeeman-Hilbert echoes retain observed density matrices,
% spectral residuals and absolute forward-error bounds in aux. The response
% remains time-by-detector, or column-vectorized observed density by time.
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

% Reuse one Hermitian spectrum and retain its forward-error diagnostics
if ~projected&&strcmp(spin_system.formalism,'zeeman-hilb')&&...
   ~any(R(:))&&~any(K(:))&&~isfield(parameters,'pump')&&...
   isfield(parameters,'cce_aux')
    rotation=expm(-1i*pi*kron((a*b'+b*a')/2,eye(db)));
    [E,aux]=generalized_echo(spin_system,parameters,H,rho0,coil,meta,rotation,density);
    return
end

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
        value=trace(rho_b);
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

% Propagate a generalized Hahn echo with retained observed density and bounds
function [answer,aux]=generalized_echo(system,parameters,H,rho0,coil,meta,pulse,density)

% Diagonalize the complete retained-cluster Hamiltonian once
H=(H+H')/2; [vectors,energies]=eig(full(H),'vector'); energies=real(energies);
residual=norm(H*vectors-vectors.*energies.','fro');
orthogonality=norm(vectors'*vectors-eye(size(H)),'fro');
dc=system.observed_dimension; db=system.bath_dimension; count=numel(parameters.times);
detectors=cellfun(@(x)x(1:db:end,1:db:end),coil,'UniformOutput',false);
answer=complex(zeros(count,numel(meta.normalization_divisor)));
if density, answer=complex(zeros(count,dc^2)); end
aux.density=complex(zeros(dc,dc,count)); bounds=zeros(size(answer));

% Spectral and detector norms are constant over the requested time mesh
operations=20*size(H,1)^3; gamma=operations*eps/(1-operations*eps);
hnorm=norm(H,'fro'); hbound=parameters.cce_aux.roundoff_bound;
detector_norms=cellfun(@(x)norm(x,2),detectors);
if ~isfield(parameters,'coil')&&~isfield(parameters,'detection'), detector_norms=1; end
divisors=meta.normalization_divisor;

% Compose both free intervals and the physical pi rotation in time order
for n=1:count
    time=parameters.times(n); propagator=eye(size(H));
    if time>0
        free=(vectors.*exp(-1i*energies*(time/2)).')*vectors';
        propagator=free*propagator; propagator=pulse*propagator;
        propagator=free*propagator;
    end
    rho=propagator*rho0*propagator'; observed=complex(zeros(dc));
    for k=1:db
        indices=k:db:dc*db; observed=observed+rho(indices,indices);
    end
    aux.density(:,:,n)=observed;
    if density
        answer(n,:)=observed(:).';
    else
        if ~isfield(parameters,'coil')&&~isfield(parameters,'detection')
            answer(n)=exp(-1i*meta.receiver_phase)*(meta.transition.alpha'*...
                      observed*meta.transition.beta)/divisors;
        else
            for k=1:numel(detectors), answer(n,k)=trace(detectors{k}'*observed)/divisors(k); end
        end
    end

    % Bound spectral residuals, tensor assembly and floating-point products
    bound=2*time*(residual+hbound)+...
          8*orthogonality+gamma*(1+time*hnorm);
    if density
        bounds(n,:)=bound;
    else
        for k=1:numel(meta.normalization_divisor)
            bounds(n,k)=bound*detector_norms(k)/abs(divisors(k));
        end
    end
end

% Keep generator and propagation evidence alongside the natural response
aux.abs_error_bound=bounds; aux.eigen_residual=residual;
aux.orthogonality_residual=orthogonality;
aux.hamiltonian_bound=parameters.cce_aux.roundoff_bound;
aux.backend=parameters.cce_aux.backend;

end


