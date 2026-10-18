% Chemical reaction network reporting. Syntax:
%
%                  summary_kinetics(spin_system)
%
% Parameters:
%
%    spin_system - Spinach data structure
%
% Outputs:
%
%    This function prints the reaction network and atom transport.
%
% i.kuprov@soton.ac.uk

function summary_kinetics(spin_system)

% Check consistency
grumble(spin_system);

% Report the reaction network and atom transport
summary_chemistry(spin_system);
for n=1:numel(spin_system.chem.reactions)
    reaction=spin_system.chem.reactions{n};
    source_spins=cellfun(@(x)x(:)',spin_system.chem.parts(reaction.reactants),'UniformOutput',false);
    source_spins=[source_spins{:}];
    report(spin_system,['reaction ' num2str(n) ': matched spins ' ...
           mat2str(reaction.matching(:,1)') ', traced spins ' ...
           mat2str(setdiff(source_spins,reaction.matching(:,1))) ...
           ', closure ' reaction.closure]);
end

end

% Consistency enforcement
function grumble(spin_system)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
end

