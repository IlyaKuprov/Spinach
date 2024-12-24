% A shorthand for inflating cell arrays of polyadics; this
% inflates every component of the cell array. Syntax:
%
%                        A=inflate(A)
%
% Parameters:
%
%    A   -  a cell array of polyadics
%
% Outputs
%
%    A   -  a cell array of matrices
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=cell/inflate.m>

function A=inflate(A)

% Inflate element by element
for n=1:numel(A)
    A{n}=inflate(A{n});
end

end

% Redemption, but not repentance.
%
% Michael Krug

