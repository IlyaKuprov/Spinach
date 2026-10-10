% Extract concentrations from a concentration-weighted state. Syntax:
%
%                  concs=chem_concs(spin_system,eta)
%
% Parameters:
%
%    spin_system - Spinach object with a compiled density-matrix basis
%
%    eta         - column vector in space-times-spin order; its length
%                  is an integer multiple of bas.offsets(end); in
%                  zeeman-hilb, a block-diagonal density matrix instead
%
% Outputs:
%
%    concs       - nvoxels-by-nsubst array of concentrations, in the
%                  concentration units used for inter.chem.concs
%
% No concentration division, clipping, or renormalisation is performed.
%
% ilya.kuprov@weizmann.ac.il

function concs=chem_concs(spin_system,eta)

% Check consistency
grumble(spin_system,eta);

% Read matrix traces directly in the Hilbert representation
if strcmp(spin_system.bas.formalism,'zeeman-hilb')
    concs=zeros(1,spin_system.bas.nsubst);
    for n=1:spin_system.bas.nsubst
        idx=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
        concs(n)=trace(eta(idx,idx));
    end
    return
end

% Read one trace functional per substance in each spatial block
eta=reshape(eta,spin_system.bas.offsets(end),[]);
if strcmp(spin_system.bas.formalism,'sphten-liouv')
    concs=eta(spin_system.bas.offsets(1:end-1)+1,:).';
else
    concs=zeros(size(eta,2),spin_system.bas.nsubst);
    for n=1:spin_system.bas.nsubst
        idx=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
        dim=prod(spin_system.comp.mults(spin_system.chem.parts{n}));
        concs(:,n)=(reshape(speye(dim),1,[])*eta(idx,:)).';
    end
end

end

% Consistency enforcement
function grumble(spin_system,eta)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv','zeeman-hilb'})
    error('Spinach:chem_concs:formalism','concentrations are not defined in zeeman-wavef formalism.');
end
if strcmp(spin_system.bas.formalism,'zeeman-hilb')
    if ~isnumeric(eta)||~isequal(size(eta),[spin_system.bas.offsets(end) spin_system.bas.offsets(end)])||...
       any(~isfinite(eta),'all')
        error('Spinach:chem_concs:state','eta must be a finite matrix of the compiled Hilbert dimensions.');
    end
    for n=1:spin_system.bas.nsubst
        idx=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
        if nnz(eta(idx,:))~=nnz(eta(idx,idx))
            error('Spinach:chem_concs:crossSubstance','eta must not contain inter-substance coherences.');
        end
    end
elseif ~isnumeric(eta)||~iscolumn(eta)||isempty(eta)||...
       mod(numel(eta),spin_system.bas.offsets(end))~=0||any(~isfinite(eta))
    error('Spinach:chem_concs:state','eta must be a finite column with an integer number of spin blocks.');
end
end


