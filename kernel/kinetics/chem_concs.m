% Extract concentrations from a concentration-weighted state. Syntax:
%
%                  concs=chem_concs(spin_system,eta)
%
% Parameters:
%
%    spin_system - Spinach object with a compiled sphten-liouv basis
%
%    eta         - column vector in space-times-spin order; its length
%                  is an integer multiple of bas.offsets(end)
%
% Outputs:
%
%    concs       - nvoxels-by-nsubst array of unit coordinates, in the
%                  concentration units used for inter.chem.concs
%
% No concentration division, clipping, or renormalisation is performed.
%
% ilya.kuprov@weizmann.ac.il

function concs=chem_concs(spin_system,eta)

% Check consistency
grumble(spin_system,eta);

% Read one unit coordinate per substance in each spatial block
eta=reshape(eta,spin_system.bas.offsets(end),[]);
concs=eta(spin_system.bas.offsets(1:end-1)+1,:).';

end

% Consistency enforcement
function grumble(spin_system,eta)
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('Spinach:chem_concs:formalism','chem_concs currently requires sphten-liouv formalism.');
end
if ~isnumeric(eta)||~iscolumn(eta)||isempty(eta)||...
   mod(numel(eta),spin_system.bas.offsets(end))~=0||any(~isfinite(eta))
    error('Spinach:chem_concs:state','eta must be a finite column with an integer number of spin blocks.');
end
end


