% Fidelity, gradient, and Hessian contributions of one block of ensem-
% ble cases, evaluated on the parallel pool worker that holds the drift
% generators of that block. This function is called by ensemble.m in-
% side its spmd block; the per-case physics (phase cycle, offsets, po-
% wer level, waveform distortions, GRAPE) is applied here. Syntax:
%
%     [traj,fid,grad,hess]=ens_block(spin_system,drifts,control,...
%                                    block,waveform,n_outputs)
%
% Parameters:
%
%   spin_system  - frozen problem published by optimcon.m, with
%                  the drift generators removed
%
%   drifts       - cell array of drift generators, populated at
%                  the indices that the cases of this block use
%
%   control      - live client-side control structure
%
%   block        - index of the case block, into the cell array
%                  spin_system.control.worker_cases
%
%   waveform     - control coefficients for each control opera-
%                  tor, [ncontrols x nsteps], rad/s
%
%   n_outputs    - number of outputs requested from ensemble.m,
%                  2 for the fidelity, 3 for the gradient, 4 for
%                  the Hessian
%
% Outputs:
%
%   traj         - cell array of trajectory structures, one per
%                  case of the block in block order; when the
%                  control.traj_opts contains 'average', one
%                  structure holding the sum over the block, or
%                  an empty cell for an empty block
%
%   fid          - [1 x n_block] array of case fidelities
%
%   grad         - sum of the case gradients over the block, a
%                  [ncontrols*nsteps x 1] column, empty unless
%                  n_outputs>2
%
%   hess         - sum of the case Hessians over the block, a
%                  [(ncontrols*nsteps)^2 x 1] column, empty un-
%                  less n_outputs>3
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ens_block.m>

function [traj,fid,grad,hess]=ens_block(spin_system,drifts,control,block,waveform,n_outputs)

% Check consistency
grumble(spin_system,drifts,control,block,waveform,n_outputs);

% Graft live client data over the frozen worker copy, keep what only the worker holds
frozen=spin_system.control; spin_system.control=control;
missing=setdiff(fieldnames(frozen),fieldnames(control));
for k=1:numel(missing)
    spin_system.control.(missing{k})=frozen.(missing{k});
end

% GRAPE function for the formalism
switch spin_system.bas.formalism
    case {'sphten-liouv','zeeman-liouv','zeeman-wavef'}
        grape=@grape_liouv;
    case 'zeeman-hilb'
        grape=@grape_hilb;
    otherwise
        error('unrecognised formalism specification.');
end

% Cases of this block, waveform dimensions, and offset ensemble size
my_cases=frozen.worker_cases{block}; n_mine=numel(my_cases); catalog=control.catalog;
ncont=size(waveform,1); nsteps=size(waveform,2); off_ens_sizes=cellfun(@numel,control.offsets);

% Preallocate block outputs, derivative buffers only when requested
traj=cell(n_mine,1); fid=zeros(1,n_mine); grad=[]; hess=[];
if n_outputs>2, grad=zeros(ncont*nsteps,1); end
if n_outputs>3, hess=zeros((ncont*nsteps)^2,1); end

% Loop over the cases of the block
for m=1:n_mine

    % Extract ensemble indices
    n=my_cases(m); n_rho=catalog(n,1); n_sys=catalog(n,2);
    n_pwr=catalog(n,3); n_off=catalog(n,4);
    n_phi=catalog(n,5); n_dis=catalog(n,6);

    % Get initial and target states, drift, and waveform
    rho_init=control.rho_init{n_rho}; rho_targ=control.rho_targ{n_rho};
    L=drifts{n_sys}; local_waveform=waveform;

    % Phase cycle: a rotation of each control channel and phases on the states
    R=eye(ncont);
    if ~isempty(control.phase_cycle)
        phi=control.phase_cycle(n_phi,:);
        rho_init=exp(1i*phi(1))*rho_init; rho_targ=exp(1i*phi(end))*rho_targ;
        R=kron(diag(cos(phi(2:end-1))),eye(2))+kron(diag(sin(phi(2:end-1))),[0 -1; 1 0]);
        local_waveform=R*local_waveform;
    end

    % Add offset terms, first channel index fastest (user specifies offsets in Hz)
    if ~isempty(off_ens_sizes)
        off_idx=cell(1,numel(off_ens_sizes)); [off_idx{:}]=ind2sub([off_ens_sizes 1],n_off);
        for k=1:numel(off_ens_sizes)
            L=L+sparse(2*pi*control.offsets{k}(off_idx{k})*spin_system.control.off_ops{k});
        end
    end

    % Move the waveform into physical units
    power_lvl=control.pwr_levels(n_pwr); local_waveform=power_lvl*local_waveform;

    % Apply waveform distortions, with their Jacobian when derivatives are needed
    if n_outputs>2, J=speye(numel(local_waveform)); end
    for k=1:size(control.distortion,2)
        if n_outputs>2
            [local_waveform,stage_jacobian]=control.distortion{n_dis,k}(local_waveform);
            J=stage_jacobian*J;
        else
            local_waveform=control.distortion{n_dis,k}(local_waveform);
        end
    end

    % Fidelity, trajectory, and derivatives
    outputs=cell(1,n_outputs);
    [outputs{:}]=grape(spin_system,L,spin_system.control.operators,local_waveform,rho_init,rho_targ,control.fidelity);
    traj{m}=outputs{1}; fid(m)=outputs{2};

    % Gradient through the Jacobian, the phase cycle, and the power level
    if n_outputs>2
        grad_n=R'*reshape(J'*outputs{3}(:),ncont,nsteps);
        grad=grad+power_lvl*grad_n(:);
    end

    % Hessian through the phase cycle and the power level
    if n_outputs>3
        K=kron(speye(nsteps),sparse(R')); hess_n=K*outputs{4}*K';
        hess=hess+power_lvl^2*hess_n(:);
    end

end

% Collapse a non-empty block into one trajectory sum when only the average is needed
if ismember('average',control.traj_opts)&&(n_mine>0)
    traj_sum=traj{1}.forward;
    for m=2:n_mine
        traj_sum=traj_sum+traj{m}.forward;
    end
    traj={struct('forward',{traj_sum})};
end

end

% Consistency enforcement
function grumble(spin_system,drifts,control,block,waveform,n_outputs)
if (~isfield(spin_system,'control'))||(~isfield(spin_system.control,'worker_cases'))
    error('spin_system must be the frozen problem published by optimcon().');
end
if ~iscell(drifts)
    error('drifts must be a cell array of drift generators.');
end
if ~isstruct(control)
    error('control must be the live control structure from ensemble().');
end
if (~isnumeric(block))||(~isscalar(block))||(mod(block,1)~=0)||...
   (block<1)||(block>numel(spin_system.control.worker_cases))
    error('block must be a positive integer not exceeding the number of case blocks.');
end
if (~isnumeric(waveform))||(~isreal(waveform))
    error('waveform must be an array of real numbers.');
end
nsteps=spin_system.control.pulse_nsteps+strcmp(spin_system.control.integrator,'trapezium');
if ~isequal(size(waveform),[spin_system.control.ncontrols nsteps])
    error('waveform must be a [ncontrols x nsteps] array, nsteps+1 columns for the trapezium integrator.');
end
if (~isnumeric(n_outputs))||(~isscalar(n_outputs))||(~ismember(n_outputs,[2 3 4]))
    error('n_outputs must be 2, 3, or 4.');
end
end

% The English climate, to be fair, is a lot better than
% England's reputation for it: it is the English who keep
% complaining about it, and everyone else believes them.
%
% Ilya Kuprov

