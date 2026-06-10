% Prints permutation-symmetry summary for a Spinach system. Syntax:
%
%                 summary_symmetry(spin_system,header)
%
% Parameters:
%
%    spin_system  - Spinach spin system description object
%
%    header       - a string of text to precede the summary
%
% Outputs:
%
%    this function prints to the console or to the user-specified
%    output via report.m function
%
% ilya.kuprov@weizmann.ac.il

function summary_symmetry(spin_system,header)

% Check consistency
grumble(spin_system,header);

% Print the symmetry table
report(spin_system,header);
report(spin_system,'=====================');
report(spin_system,' Group    Spins      ');
report(spin_system,'---------------------');
for n=1:length(spin_system.comp.sym_spins)
    report(spin_system,['  ' spin_system.comp.sym_group{n} '     ' num2str(spin_system.comp.sym_spins{n})]);
end
report(spin_system,'=====================');

end

% Consistency enforcement
function grumble(spin_system,header)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
if ~ischar(header)
    error('header must be a character string.');
end
end

