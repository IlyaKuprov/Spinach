% Returns the identity on the compiled direct-sum basis. Its dimension
% is bas.offsets(end), the sum of the substance dimensions. Syntax:
%
%                      A=unit_oper(spin_system)
%
% Parameters:
%
%    spin_system  - Spinach data object containing basis 
%                   information (call basis.m first)
%
% Outputs:
%
%    A            - a sparse unit matrix of appropriate
%                   dimension 
% 
% ilya.kuprov@weizmann.ac.il
% d.savostyanov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=unit_oper.m>

function A=unit_oper(spin_system)

% Check consistency
grumble(spin_system);

% Unit matrix on the direct sum
A=speye(spin_system.bas.offsets(end));

end

% Consistency enforcement
function grumble(spin_system)
if (~isfield(spin_system,'bas'))||(~isfield(spin_system.bas,'formalism'))
    error('the spin_system object does not contain the required information.');
end
end

% The substance of this book, as it is expressed in the editor's preface, is
% that to measure "right" by the false philosophy of the Hebrew prophets and
% "weepful" Messiahs is madness. Right is not the offspring of doctrine, but
% of power. All laws, commandments, or doctrines as to not doing to another
% what you do not wish done to you, have no inherent authority whatever, but
% receive it only from the club, the gallows and the sword. A man truly free
% is under no obligation to obey any injunction, human or divine. [...] Men
% should not be bound by moral rules invented by their foes.
%
% Leo Tolstoy, about Ragnar Redbeard's "Might is Right"

