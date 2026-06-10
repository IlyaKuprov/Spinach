% Prints periodic boundary condition vector summary for a Spinach system. Syntax:
%
%                 summary_pbc(spin_system,header)
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

function summary_pbc(spin_system,header)

% Check consistency
grumble(spin_system,header);

% Print the vector table
report(spin_system,header);
report(spin_system,'===============================');
report(spin_system,'     X         Y         Z     ');
report(spin_system,'-------------------------------');
for n=1:numel(spin_system.inter.pbc)
    report(spin_system,['   ' pad(num2str(spin_system.inter.pbc{n}(1),'%+5.3f   '),10)...
                              pad(num2str(spin_system.inter.pbc{n}(2),'%+5.3f   '),10)...
                              pad(num2str(spin_system.inter.pbc{n}(3),'%+5.3f   '),10)]);
end
report(spin_system,'===============================');

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

