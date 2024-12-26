% Dot product of TT representations of matrices. Syntax:
%
%                        c=dot(a,b)
%
% Parameters:
%
%    a,b  - tensor train objects representing numerical
%           arrays of consistent dimensions and having
%           the same internal topology
%
% Outputs:
%
%    c    - inner product of a and b
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/dot.m>

function c=dot(a,b)

% Check consistency
grumble(a,b);

% Compute the product
c=ctranspose(a)*b;

end

% Consistency enforcement
function grumble(a,b)
if (~isa(a,'ttclass'))||(~isa(b,'ttclass'))
    error('both inputs should be tensor trains.')
end
if (size(a.cores,1)~=size(b.cores,1))||...
   (~all(all(sizes(a)==sizes(b))))
    error('tensor train topologies must be the same.')
end
end

% Any product that needs a manual to work is broken.
%
% Elon Musk

