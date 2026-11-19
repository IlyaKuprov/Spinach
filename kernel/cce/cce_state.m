% Construct a state in the declared retained-cluster representation. Syntax:
%
%                  rho=cce_state(cluster_system,states,spins)
%
% Character descriptors and cell rows of character descriptors follow
% cce_operator, with conventional unnormalised spin-density operators.
% Explicit physical states use a cell row of local density matrices or
% population vectors, and a matching numeric physical-ID row. Every retained
% spin must occur exactly once: no unspecified population is invented.
% Population vectors follow descending magnetic quantum number. Matrices
% must be positive, Hermitian, and trace one. Their tensor-product order is
% always cluster_system.physical_ids, not the input descriptor order.
% Output is a sparse Hilbert density matrix or its column-major Liouville
% vector. No full-bath state is built, and physical states are not rescaled.
%
% Parameters:
%
%   cluster_system - retained subsystem returned by cce_subsystem
%   states         - operator descriptors or a cell row of local physical
%                    density matrices/population vectors
%   spins          - matching operator selections or a physical-ID row
%
% Outputs:
%
%   rho - sparse retained density matrix or column-major Liouville vector
%
% ilya.kuprov@weizmann.ac.il

function rho=cce_state(cluster_system,states,spins)

% Validate descriptors before changing the local representation
grumble(cluster_system,states,spins);
formalism=cluster_system.formalism;
cluster_system.formalism='zeeman-hilb';
if ischar(states)||(iscell(states)&&all(cellfun(@ischar,states)))

    % Conventional spin-density operators use the common operator algebra
    rho=cce_operator(cluster_system,states,spins);
else

    % Every physical population is supplied explicitly in local dimension
    ids=cluster_system.physical_ids;
    rho=sparse(1);
    [~,order]=ismember(ids,spins);
    numerics=cluster_system.model.numerics;
    for n=1:numel(ids)
        d=cluster_system.multiplicities(n);
        local=states{order(n)};
        if isvector(local)&&numel(local)==d
            local=diag(local);
        end
        if 16*(size(rho,1)*d)^2>numerics.max_memory_bytes
            error('cce_state:memory','Requested density exceeds the memory budget.');
        end
        rho=kron(rho,sparse(local));
    end

end

% Use the same column-major vectorisation as stock Spinach states
if strcmp(formalism,'zeeman-liouv')
    rho=rho(:);
end

end

% Keep operator descriptors separate from explicit physical populations
function grumble(system,states,spins)
if ~isstruct(system)||~isscalar(system)||~isfield(system,'kind')||...
   ~strcmp(system.kind,'cce_subsystem')||...
   ~ismember(system.formalism,{'zeeman-hilb','zeeman-liouv'})
    error('cce_state:system','Use a supported cce_subsystem view.');
end
if ischar(states)||(iscell(states)&&all(cellfun(@ischar,states)))
    if ischar(states)
        states={states};
        spins={spins};
    end
    if ~isrow(states)||isempty(states)||~all(cellfun(@isrow,states))||...
       ~iscell(spins)||~isequal(size(states),size(spins))
        error('cce_state:descriptor','Supply a label or matching cell rows of labels and spin selections.');
    end
    for n=1:numel(states)
        if ~ismember(states{n},{'E','Lx','Ly','Lz','L+','L-'})
            error('cce_state:descriptor','Unsupported operator descriptor.');
        end
        ids=spins{n};
        if ischar(ids)
            if strcmp(ids,'all'), ids=system.physical_ids;
            else
                ids=system.physical_ids(strcmp(system.isotopes,ids));
            end
        end
        if ~isnumeric(ids)||~isreal(ids)||isempty(ids)||~isrow(ids)||...
           numel(unique(ids))~=numel(ids)||any(~ismember(ids,system.physical_ids))
            error('cce_state:descriptor','Select distinct retained physical IDs.');
        end
    end
    return
end
if ~iscell(states)||~isrow(states)||~isnumeric(spins)||~isrow(spins)||...
   numel(states)~=numel(spins)||~isequal(sort(spins),sort(system.physical_ids))
    error('cce_state:descriptor','Explicit states must specify every retained physical spin once.');
end
[~,order]=ismember(system.physical_ids,spins);
numerics=system.model.numerics;
for n=1:numel(system.physical_ids)
    d=system.multiplicities(n);
    local=states{order(n)};
    if isvector(local)&&numel(local)==d
        local=diag(local);
    end
    tolerance=numerics.atol+numerics.rtol*norm(local,'fro');
    if ~isnumeric(local)||~isequal(size(local),[d d])||...
       any(~isfinite(local),'all')||norm(local-local','fro')>tolerance||...
       abs(trace(local)-1)>tolerance||...
       min(real(eig(full((local+local')/2))))<-tolerance
        error('cce_state:density','Supply positive, Hermitian, trace-one local density matrices.');
    end
end
end


