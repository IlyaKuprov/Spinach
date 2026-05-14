% Nucleus record builder for diamond defects. Syntax:
%
%          nucleus=diamond_nuc(iso,A,Q)
%
% This internal helper builds one nucleus record for diamond-defect
% constructors.
%
% Parameters:
%
%    iso      - Spinach isotope label
%
%    A        - hyperfine tensor in Hz
%
%    Q        - quadrupole tensor in Hz, or []
%
% Outputs:
%
%    nucleus  - isotope, hyperfine, and quadrupole record
%
% <https://spindynamics.org/wiki/index.php?title=diamond_nuc.m>

function nucleus=diamond_nuc(iso,A,Q)

% Set default quadrupole tensor
if nargin<3
    Q=[];
end

% Check consistency
grumble(iso,A,Q);

% Build the nucleus record
nucleus=struct('iso',iso,'A',A,'Q',Q);

end

% Consistency enforcement
function grumble(iso,A,Q)
if ~ischar(iso)
    error('iso must be a character string.');
end
if(~isnumeric(A)||~isreal(A)||~isequal(size(A),[3 3]))
    error('A must be a real 3x3 matrix.');
end
if(~isempty(Q)&&(~isnumeric(Q)||~isreal(Q)||~isequal(size(Q),[3 3])))
    error('Q must be empty or a real 3x3 matrix.');
end
end

% A nucleus record should carry only what the Hamiltonian needs.

