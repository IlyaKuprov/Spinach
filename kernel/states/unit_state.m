% Returns the unit state vector or matrix in the current formalism
% and basis, weighted by each substance concentration. Syntax:
%
%                     rho=unit_state(spin_system)
%
% Parameters:
%
%    spin_system  - Spinach data object containing basis 
%                   information (call basis.m first)
%
% Outputs:
%
%    rho          - vector or matrix representation of
%                   the unit state 
%
% Note: Zeeman blocks retain the stock geometric identity normalisation.
%       Use equilibrium for trace-one density matrices before weighting.
%
% ilya.kuprov@weizmann.ac.il
% d.savostyanov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=unit_state.m>

function rho=unit_state(spin_system)

% Check consistency
grumble(spin_system);

% Decide how to proceed
switch spin_system.bas.formalism
    
    case 'sphten-liouv'
        
        % Concentration at each T(0,0) coordinate
        rho=sparse(spin_system.bas.offsets(1:end-1)+1,1,spin_system.chem.concs,...
                   spin_system.bas.offsets(end),1);
        
    case 'zeeman-liouv'
        
        % Stack locally normalised stretched identities
        blocks=cell(spin_system.bas.nsubst,1);
        for n=1:spin_system.bas.nsubst
            block=speye(prod(spin_system.comp.mults(spin_system.chem.parts{n})));
            block=block(:); blocks{n}=spin_system.chem.concs(n)*block/norm(block,2);
        end
        rho=vertcat(blocks{:});
        
    case 'zeeman-hilb'
        
        % Place weighted geometric identities on the Hilbert diagonal
        blocks=cell(spin_system.bas.nsubst,1);
        for n=1:spin_system.bas.nsubst
            blocks{n}=spin_system.chem.concs(n)*...
                      speye(prod(spin_system.comp.mults(spin_system.chem.parts{n})));
        end
        rho=blkdiag(blocks{:});
        
    otherwise
        
        % Complain and bomb out
        error('unknown formalism specification.');
        
end

end

% Consistency enforcement
function grumble(spin_system)
if (~isfield(spin_system,'bas'))||(~isfield(spin_system.bas,'formalism'))
    error('the spin_system object does not contain the required information.');
end
if strcmp(spin_system.bas.formalism,'zeeman-wavef')
    error('Spinach:unit_state:wavefunction',...
          'concentration-weighted unit states are not supported in zeeman-wavef formalism.');
end
end

% There used to be a simple story about Russian literature, that we
% thought the good writers were the ones who opposed the regime. Once
% we don't have that story about Russia as a competitor, or an enemy,
% it was much less clear to us what we should be interested in.
%
% Edwin Frank, the editor of NYRB Classics

