% Reconstructs constructor descriptions of the cores of a polyadic.
% Syntax:
%
%                         cores=core_specs(p)
%
% Parameters:
%
%    p     - a polyadic object
%
% Outputs:
%
%    cores - nested core cells accepted by the polyadic constructor;
%            implicit actions include their adjoints and dimensions
%
% ilya.kuprov@weizmann.ac.il

function cores=core_specs(p)

% Check consistency
grumble(p);

% Pair implicit actions with their construction metadata
cores=p.cores;
for n=1:numel(cores)
    for k=1:numel(cores{n})
        if isa(cores{n}{k},'function_handle')
            cores{n}{k}=struct('action',cores{n}{k},...
                              'adjoint',p.core_adj{n}{k},...
                              'dims',p.core_dims{n}{k});
        end
    end
end

end

% Consistency enforcement
function grumble(p)
if ~isa(p,'polyadic')
    error('p must be polyadic.');
end
end


