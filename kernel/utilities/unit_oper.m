% Returns a unit operator in the current formalism and basis. The
% operator has dimension equal to the basis size in sphten-liouv
% formalism, the dimension equal to the product of all spin multi-
% plicities in zeeman-hilb formalism, and the dimension of square
% of the product of all spin multiplicities in zeeman-liouv forma-
% lism. Syntax:
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
% i.kuprov@soton.ac.uk
% d.savostyanov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=unit_oper.m>

function A=unit_oper(spin_system)

% Check consistency
grumble(spin_system);

% Decide how to proceed
switch spin_system.bas.formalism
    
    case 'sphten-liouv'
        
        % Unit matrix
        A=speye(size(spin_system.bas.basis,1));
        
    case 'zeeman-liouv'
        
        % Unit matrix
        A=speye(prod(spin_system.comp.mults.^2));
        
    case 'zeeman-hilb'
        
        % Unit matrix
        A=speye(prod(spin_system.comp.mults));
        
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

