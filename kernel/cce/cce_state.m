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
% must be positive, Hermitian and trace one. Their tensor-product order is
% always cluster_system.physical_ids, not the input descriptor order.
% Output is a sparse Hilbert density matrix or its column-major Liouville
% vector. No full-bath state is built, and physical states are not rescaled.
%
% ilya.kuprov@weizmann.ac.il

function rho=cce_state(cluster_system,states,spins)

% Validate descriptors before changing the local representation
grumble(cluster_system,states,spins);
formalism=cluster_system.formalism; cluster_system.formalism='zeeman-hilb';
if ischar(states)||(iscell(states)&&all(cellfun(@ischar,states)))

    % Conventional spin-density operators use the common operator algebra
    rho=cce_operator(cluster_system,states,spins);
else

    % Every physical population is supplied explicitly in local dimension
    ids=cluster_system.physical_ids; rho=sparse(1);
    [~,order]=ismember(ids,spins); numerics=cluster_system.model.numerics;
    for n=1:numel(ids)
        d=cluster_system.multiplicities(n); local=states{order(n)};
        if isvector(local)&&numel(local)==d, local=diag(local); end
        tolerance=numerics.atol+numerics.rtol*norm(local,'fro');
        if ~isnumeric(local)||~isequal(size(local),[d d])||...
           any(~isfinite(local),'all')||norm(local-local','fro')>tolerance||...
           abs(trace(local)-1)>tolerance||...
           min(real(eig(full((local+local')/2))))<-tolerance
            error('cce_state:density','Supply positive, Hermitian, trace-one local density matrices.');
        end
        if 16*(size(rho,1)*d)^2>numerics.max_memory_bytes
            error('cce_state:memory','Requested density exceeds the memory budget.');
        end
        rho=kron(rho,sparse(local));
    end
end

% Use the same column-major vectorisation as stock Spinach states
if strcmp(formalism,'zeeman-liouv'), rho=rho(:); end

end

% Keep operator descriptors separate from explicit physical populations
function grumble(system,states,spins)
if ~isstruct(system)||~isscalar(system)||~isfield(system,'kind')||...
   ~strcmp(system.kind,'cce_subsystem')||...
   ~ismember(system.formalism,{'zeeman-hilb','zeeman-liouv'})
    error('cce_state:system','Use a supported cce_subsystem view.');
end
if ischar(states)||(iscell(states)&&all(cellfun(@ischar,states))), return; end
if ~iscell(states)||~isrow(states)||~isnumeric(spins)||~isrow(spins)||...
   numel(states)~=numel(spins)||~isequal(sort(spins),sort(system.physical_ids))
    error('cce_state:descriptor','Explicit states must specify every retained physical spin once.');
end
end


