% Prints chemical subsystem and reaction summary for a Spinach system. Syntax:
%
%                 summary_chemistry(spin_system)
%
% Parameters:
%
%    spin_system  - Spinach spin system description object
%
% Outputs:
%
%    this function prints to the console or to the user-specified
%    output via report.m function
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=summary_chemistry.m>

function summary_chemistry(spin_system)

% Check consistency
grumble(spin_system);

% Report substance membership and initial concentration
for n=1:numel(spin_system.chem.parts)
    report(spin_system,['chemical subsystem ' num2str(n) ': spins ['...
                       num2str(spin_system.chem.parts{n}(:)') '], concentration '...
                       num2str(spin_system.chem.concs(n))]);
end

% Report explicit reaction records and their closures
for n=1:numel(spin_system.chem.reactions)
    reaction=spin_system.chem.reactions{n};
    if isa(reaction.rate,'function_handle')
        rate_text=func2str(reaction.rate);
    else
        rate_text=num2str(reaction.rate);
    end
    report(spin_system,['reaction ' num2str(n) ': [' num2str(reaction.reactants)...
                       '] -> [' num2str(reaction.products) '], rate ' rate_text...
                       ', closure ' reaction.closure]);
    report(spin_system,['matched spin pairs: ' mat2str(reaction.matching)]);
    if isfield(reaction,'selector')
        if ischar(reaction.selector{1})
            report(spin_system,['selector: ' reaction.selector{1}]);
        else
            report(spin_system,'selector: user product superoperator pair');
        end
    end
end

end

% Consistency enforcement
function grumble(spin_system)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
end

% Always code as if the guy who ends up maintaining 
% your code will be a violent psychopath who knows 
% where you live.
%
% Martin Golding

