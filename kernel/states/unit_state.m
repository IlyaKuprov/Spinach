% Returns the unit state vector or matrix in the current formalism
% and basis. Syntax:
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
% Note: multi-substance Zeeman unit states are not yet supported.
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
        
        % Unit population of T(0,0) state
        rho=sparse(spin_system.bas.offsets(1:end-1)+1,1,1,...
                   spin_system.bas.offsets(end),1);
        
    case 'zeeman-liouv'
        
        % Normalized stretched unit matrix
        rho=speye(prod(spin_system.comp.mults));
        rho=rho(:); rho=rho/norm(rho,2);
        
    case 'zeeman-hilb'
        
        % Sparse unit matrix
        rho=speye(prod(spin_system.comp.mults));
        
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
if ismember(spin_system.bas.formalism,{'zeeman-liouv','zeeman-hilb'})&&...
   (spin_system.bas.nsubst>1)
    error('Spinach:unit_state:segmentedZeeman',...
          'multi-substance Zeeman unit states are not yet supported.');
end
end

% There used to be a simple story about Russian literature, that we
% thought the good writers were the ones who opposed the regime. Once
% we don't have that story about Russia as a competitor, or an enemy,
% it was much less clear to us what we should be interested in.
%
% Edwin Frank, the editor of NYRB Classics

