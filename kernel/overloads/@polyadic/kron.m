% Kronecker product function for polyadics. Syntax:
%
%                      c=kron(a,b)
%
% Parameters:
%
%    a,b - polyadic or numeric objects
%
% Outputs:
%
%    c - polyadic object
%
% This operation bundles the inputs into a nested polyadic object.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/kron.m>

function c=kron(a,b)

% Check consistency
grumble(a,b);

% Put the new term inside the polyadic structure
if isa(a,'polyadic')&&(~isa(b,'polyadic'))&&isempty(a.suffix)
    
    % Append B to core lists of A
    for n=1:numel(a.cores)
        a.cores{n}=[a.cores{n} {b}]; 
    end
    c=a;
    
elseif (~isa(a,'polyadic'))&&isa(b,'polyadic')&&isempty(b.prefix)
    
    % Prepend A to core lists of B
    for n=1:numel(b.cores)
        b.cores{n}=[{a} b.cores{n}];
    end
    c=b;
    
else
    
    % Make a nested polyadic
    c=polyadic({{a,b}});
    
end

% Simplify
c=simplify(c);

end

% Consistency enforcement
function grumble(a,b)
if (~isa(a,'polyadic'))&&(~isnumeric(a))
    error('operands must be either matrices or polyadics.');
end
if (~isa(b,'polyadic'))&&(~isnumeric(b))
    error('operands must be either matrices or polyadics.');
end
end

% "Order of authorship was determined by proximity to 
%  tenure decisions."
%
% Acknowledgement in 
% http://dx.doi.org/10.1046/j.1365-294x.1998.00309.x

