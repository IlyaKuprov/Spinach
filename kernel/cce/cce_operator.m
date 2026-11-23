% Construct canonical spin operators on a retained cluster. Syntax:
%
%                 A=cce_operator(cluster_system,operators,spins)
%
% cluster_system is cce_subsystem output. operators is Lx, Ly, Lz, L+, L-
% or E. spins is a physical-ID row, isotope label, or 'all'; a row selects
% a sum. Cell rows of operators and spins select an ordered product of sums,
% as in the stock operator service. E is the dimensionless unit operator.
% Hilbert output is a sparse matrix; zeeman-liouv output is the commutation
% superoperator for column-major density vectorisation. No hidden normalisation
% or physical population is applied. Unsupported labels or absent IDs fail.
%
% Parameters:
%
%   cluster_system - retained subsystem returned by cce_subsystem
%   operators      - character label or cell row of canonical operator labels
%   spins          - physical-ID row, isotope label, 'all', or matching cell row
%
% Outputs:
%
%   A - sparse Hilbert operator or column-major commutation superoperator
%
% ilya.kuprov@weizmann.ac.il

function A=cce_operator(cluster_system,operators,spins)

% Check the view, descriptor shape, and explicit tensor-product selections
grumble(cluster_system,operators,spins);
if ischar(operators)
    operators={operators};
    spins={spins};
end

% Start from the unit operator on the retained space
d=cluster_system.cce.dimension;
A=speye(d);
labels={'E','Lx','Ly','Lz','L+','L-'};
fields={'u','x','y','z','p','m'};
for n=1:numel(operators)

    % Resolve isotope and physical-ID selections in the retained ordering
    ids=spins{n};
    if ischar(ids)
        if strcmp(ids,'all')
            ids=cluster_system.cce.ids;
        else
            ids=cluster_system.cce.ids(strcmp(cluster_system.comp.isotopes(cluster_system.cce.ids),ids));
        end
    end
    position=find(strcmp(labels,operators{n}),1);
    term=sparse(d,d);
    for id=ids

        % Embed only the requested canonical local matrix
        local=cluster_system.cce.physical_to_local(id);
        dims=cluster_system.comp.mults(cluster_system.cce.ids);
        matrix=cluster_system.cce.local_operators{local}.(fields{position});
        term=term+kron(kron(speye(prod(dims(1:local-1))),matrix),...
                       speye(prod(dims(local+1:end))));

    end
    A=A*term;

end

% Stock commutator convention for column-major density vectors
if strcmp(cluster_system.bas.formalism,'zeeman-liouv')
    A=kron(speye(d),A)-kron(A.',speye(d));
end

end

% Validate supported descriptors without constructing any operator
function grumble(system,operators,spins)
if ~isstruct(system)||~isscalar(system)||~isfield(system,'cce')||...
   ~isfield(system.cce,'sites')||...
   ~ismember(system.bas.formalism,{'zeeman-hilb','zeeman-liouv'})
    error('cce_operator:system','Use a supported cce_subsystem view.');
end
if ischar(operators)&&isrow(operators)
    operators={operators};
    spins={spins};
end
if ~iscell(operators)||~isrow(operators)||isempty(operators)||...
   ~all(cellfun(@(label)ischar(label)&&isrow(label),operators))||...
   ~iscell(spins)||~isequal(size(operators),size(spins))
    error('cce_operator:descriptor','Supply a label or matching cell rows of labels and spin selections.');
end
for n=1:numel(operators)
    ids=spins{n};
    if ischar(ids)
        if strcmp(ids,'all')
            ids=system.cce.ids;
        else
            ids=system.cce.ids(strcmp(system.comp.isotopes(system.cce.ids),ids));
            if isempty(ids)
                error('cce_operator:spins','No retained spin has this isotope.');
            end
        end
    end
    if ~isnumeric(ids)||~isreal(ids)||isempty(ids)||~isrow(ids)||...
       numel(unique(ids))~=numel(ids)||any(~ismember(ids,system.cce.ids))
        error('cce_operator:spins','Select distinct retained physical IDs.');
    end
    if ~ismember(operators{n},{'E','Lx','Ly','Lz','L+','L-'})
        error('cce_operator:label','Unsupported operator descriptor.');
    end
end
end


