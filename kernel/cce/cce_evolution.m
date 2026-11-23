% Propagate a retained cluster using conventional or projected branch generators.
% Syntax:
%
%   [answer,state_error]=cce_evolution(spin_system,L,coil,rho,timestep,nsteps,output)
%
% L is H+1i*R+1i*K for conventional generators, or the projected cell
% {left,right,dissipator}. timestep is seconds; nsteps is a nonnegative integer. 'final'
% returns the state after nsteps intervals; 'trajectory' includes t=0 in
% its columns (Liouville/projected) or pages (Hilbert); 'observable' returns
% a column per detector, including t=0, with Hilbert-Schmidt detection.
% 'total' sums trajectory states, including t=0. No reference division or
% receiver phase is applied here. Projected input/output is a column-major
% bath-coherence vector. Conventional Hilbert states are square matrices.
% The small retained-space implementation avoids a repeated basis build.
%
% Parameters:
%
%   spin_system - retained subsystem returned by cce_subsystem
%   L           - angular-frequency generator or projected cell
%   coil        - detector matrix/vector or cell of detectors
%   rho         - initial Hilbert matrix or column-major state vector
%   timestep    - finite time step in seconds
%   nsteps      - nonnegative integer number of intervals
%   output      - 'final', 'trajectory', 'observable', or 'total'
%
% Outputs:
%
%   answer - propagated states or detector values in the layout above
%
%   state_error - optional column of absolute Frobenius-norm propagation
%                 estimates for the evaluated states (one for 'final').
%                 These measure floating-point propagation only, before
%                 detection or summation; they are not rigorous bounds.
%
% ilya.kuprov@weizmann.ac.il

function [answer,state_error]=cce_evolution(spin_system,L,coil,rho,timestep,nsteps,output)

% Reject unsupported representation/output pairs before propagation
grumble(spin_system,L,coil,rho,timestep,nsteps,output);
projected=iscell(L);
hilbert=~projected&&strcmp(spin_system.bas.formalism,'zeeman-hilb');
if projected
    d=spin_system.cce.bath_dimension;
else
    d=spin_system.cce.dimension;
end
if iscell(coil)
    coils=coil;
else
    coils={coil};
end

% Keep the established eigenvalue propagation for coherent retained systems
if ~projected&&~hilbert
    L=full(L);
end
if projected
    [va,ea]=eig(full((L{1}+L{1}')/2),'vector');
    [vb,eb]=eig(full((L{2}+L{2}')/2),'vector');
    if ~isempty(L{3})
        generator=-1i*(kron(eye(d),L{1})-kron(L{2}.',eye(d)))+L{3};
    end
elseif hilbert
    [vectors,energies]=eig(full((L+L')/2),'vector');
end

% Allocate only the requested output, preserving stock trajectory orientation
switch output
    case 'observable'
        answer=complex(zeros(nsteps+1,numel(coils)));
    case 'trajectory'
        if hilbert
            answer=complex(zeros(d,d,nsteps+1));
        else
            answer=complex(zeros(d^2,nsteps+1));
        end
    case 'total'
        answer=zeros(size(rho));
    case 'final'
        answer=rho;
end
steps=0:nsteps;
if strcmp(output,'final')
    steps=nsteps;
end

% Estimate state roundoff using the CCE segment arithmetic models
if nargout>1
    elapsed=abs(timestep)*steps(:);
    if projected&&isempty(L{3})
        frequency=norm(L{1},1)+norm(L{2},1);
        roundoff=16*d^3*eps;
    elseif hilbert
        frequency=2*norm(L,1);
        roundoff=16*d^3*eps;
    else
        if projected
            frequency=norm(generator,1);
        else
            frequency=norm(L,1);
        end
        roundoff=64*d^2*eps;
    end
    state_error=roundoff*norm(rho,'fro')*((elapsed>0)+elapsed*frequency);
end

% Propagate each requested state without compounding previous time samples
for step=steps

    % Evaluate from the original state to avoid cumulative time-step drift
    time=step*timestep;
    if time==0
        state=rho;
    elseif projected
        if isempty(L{3})
            ua=va*diag(exp(-1i*ea*time))*va';
            ub=vb*diag(exp(-1i*eb*time))*vb';
            state=ua*reshape(rho,d,d)*ub';
            state=state(:);
        else
            state=expm(generator*time)*rho;
        end
    elseif hilbert
        propagator=vectors*diag(exp(-1i*energies*time))*vectors';
        state=propagator*rho*propagator';
    else
        state=expm(-1i*L*time)*rho;
    end

    % Detection is conjugate-linear in the coil, including non-Hermitian coils
    switch output
        case 'observable'
            for k=1:numel(coils)
                answer(step+1,k)=coils{k}(:)'*state(:);
            end
        case 'trajectory'
            if hilbert
                answer(:,:,step+1)=state;
            else
                answer(:,step+1)=state;
            end
        case 'total'
            answer=answer+state;
        case 'final'
            answer=state;
    end

end

end

% Validate dimensions and explicit evolution semantics
function grumble(system,L,coil,rho,timestep,nsteps,output)
if ~isstruct(system)||~isfield(system,'cce')||~isfield(system.cce,'sites')
    error('cce_evolution:system','A native retained subsystem is required.');
end
if ~isnumeric(timestep)||~isreal(timestep)||~isscalar(timestep)||~isfinite(timestep)||...
   ~isnumeric(nsteps)||~isreal(nsteps)||~isscalar(nsteps)||~isfinite(nsteps)||...
   nsteps<0||nsteps~=fix(nsteps)||~ischar(output)||...
   ~ismember(output,{'final','trajectory','observable','total'})
    error('cce_evolution:output','Use a finite time step, integer count and supported output.');
end
if ~isnumeric(rho)||any(~isfinite(rho),'all')
    error('cce_evolution:state','The initial state must be finite numeric data.');
end
if iscell(L)
    if numel(L)~=3
        error('cce_evolution:branches','Supply the two branch Hamiltonians and the dissipator.');
    end
    d=system.cce.bath_dimension;
    if ~isequal(size(rho),[d^2 1])||~isequal(size(L{1}),[d d])||...
       ~isequal(size(L{2}),[d d])||any(~isfinite(L{1}),'all')||...
       any(~isfinite(L{2}),'all')||(~isempty(L{3})&&...
       (~isequal(size(L{3}),[d^2 d^2])||any(~isfinite(L{3}),'all')))
        error('cce_evolution:branches','Projected branch dimensions must match the bath block.');
    end
    if norm(L{1}-L{1}','fro')>d*eps(max(1,norm(L{1},'fro')))||...
       norm(L{2}-L{2}','fro')>d*eps(max(1,norm(L{2},'fro')))
        error('cce_evolution:hermiticity','Hilbert branch generators must be Hermitian.');
    end
elseif strcmp(system.bas.formalism,'zeeman-hilb')
    d=system.cce.dimension;
    if ~isnumeric(L)||any(~isfinite(L),'all')||~isequal(size(L),[d d])||...
       ~isequal(size(rho),[d d])
        error('cce_evolution:hilbert','Hilbert generators and states must be square retained matrices.');
    end
    if norm(L-L','fro')>d*eps(max(1,norm(L,'fro')))
        error('cce_evolution:hermiticity','Hilbert branch generators must be Hermitian.');
    end
elseif strcmp(system.bas.formalism,'zeeman-liouv')
    d=system.cce.dimension;
    if ~isnumeric(L)||any(~isfinite(L),'all')||~isequal(size(L),[d^2 d^2])||...
       ~isequal(size(rho),[d^2 1])
        error('cce_evolution:liouville','Liouville generators act on column-major density vectors.');
    end
else
    error('cce_evolution:formalism','Only complete Zeeman representations are supported.');
end
if strcmp(output,'observable')
    if ~iscell(coil)
        coil={coil};
    end
    if isempty(coil)||~all(cellfun(@(x)isnumeric(x)&&numel(x)==d^2&&...
                                 all(isfinite(x),'all'),coil))
        error('cce_evolution:coil','Every detector must match the propagated state.');
    end
end
elements=numel(rho);
if strcmp(output,'trajectory')
    elements=elements*(nsteps+1);
end
if strcmp(output,'observable')
    elements=numel(coil)*(nsteps+1);
end
bytes=16*(12*d^2+elements);
if (iscell(L)&&~isempty(L{3}))||...
   (~iscell(L)&&strcmp(system.bas.formalism,'zeeman-liouv'))
    bytes=bytes+16*4*d^4;
end
if bytes>system.tols.cce_mem
    error('cce:memory','Propagation exceeds sys.tols.cce_mem.');
end
end


