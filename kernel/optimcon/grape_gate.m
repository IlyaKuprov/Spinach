% Evaluates squared normalised trace overlap with a target unitary.
% Syntax: [traj,fidelity,grad]=grape_gate(spin_system,drifts,ops,wave,target)
%
% The rectangle waveform is in physical angular frequency units. Drifts
% contain either one Hamiltonian or one Hamiltonian per slice. Trajectory
% output is empty because no density operator is propagated. Slice
% histories are retained only when the gradient is requested.
%
% Talos <talos@spindynamics.org>

function [traj,fidelity,grad]=grape_gate(spin_system,drifts,ops,wave,target)

% Check consistency
grumble(spin_system,drifts,ops,wave,target);

% Initialise the running time-ordered product
nsteps=size(wave,2);
dim=size(target,1);
U=eye(dim);
dt=spin_system.control.pulse_dt;
need_grad=nargout>2;

% Retain slice histories only for the backward derivative sweep
if need_grad
    prefix=cell(1,nsteps);
    H=cell(1,nsteps);
    P=cell(1,nsteps);
end

% Accumulate the pulse in chronological order
for n=1:nsteps

    % Assemble and exponentiate the physical slice
    Hslice=drifts{min(n,numel(drifts))};
    for k=1:numel(ops)
        Hslice=Hslice+wave(k,n)*ops{k};
    end
    Pslice=propagator(spin_system,Hslice,dt(n));

    % Store the product strictly before this slice for differentiation
    if need_grad
        prefix{n}=U;
        H{n}=Hslice;
        P{n}=Pslice;
    end
    U=Pslice*U;

end
overlap=trace(target'*U)/dim;
fidelity=abs(overlap)^2;
traj=struct('forward',[]);

% Differentiate the overlap using a backward suffix sweep
if need_grad

    % Move the target into the left factor for the trace
    suffix=target';
    grad=zeros(size(wave));
    for n=nsteps:-1:1

        % Contract each slice derivative with the same full overlap
        for k=1:numel(ops)
            D=dirdiff(spin_system,H{n},ops{k},dt(n),2);
            grad(k,n)=2*real(conj(overlap)*trace(suffix*D{2}*prefix{n}))/dim;
        end
        suffix=suffix*P{n};

    end

end

end

% Consistency enforcement
function grumble(spin_system,drifts,ops,wave,target)
if ~strcmp(spin_system.bas.formalism,'zeeman-hilb')||...
   ~strcmp(spin_system.control.integrator,'rectangle')
    error('gate evaluation requires rectangle Hilbert-space propagation.');
end
if ~isnumeric(wave)||~isreal(wave)||~ismatrix(wave)||...
   size(wave,1)~=numel(ops)||size(wave,2)~=numel(spin_system.control.pulse_dt)
    error('gate waveform dimensions must match the operators and time slices.');
end
if ~iscell(drifts)||~ismember(numel(drifts),[1 size(wave,2)])||...
   ~iscell(ops)||~isnumeric(target)||size(target,1)~=size(target,2)
    error('gate drift, controls, or target dimensions are invalid.');
end
end


