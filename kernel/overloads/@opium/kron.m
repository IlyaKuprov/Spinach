% Kronecker products involving an OPIUM object. Syntax:
%
%                          c=kron(a,b)
%
% Parameters:
%
%     a,b   - Kronecker operands, can be
%             matrices or opia
%
% Outputs:
%
%       c   - resulting product
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=opium/kron.m>

function c=kron(a,b)

% When both are opia
if isa(a,'opium')&&isa(b,'opium')
    
    % Return a bigger opium
    c=opium(a.dim*b.dim,a.coeff*b.coeff); return
    
end

% When A is an opium
if isa(a,'opium')
    
    % Inflate and do the kron
    c=kron(a.coeff*speye(a.dim),b); return
    
end

% When B is an opium
if isa(b,'opium')
    
    % Inflate and do the kron
    c=kron(a,b.coeff*speye(b.dim)); return
    
end

% Complain and bomb out
error('operands must be either numeric or opium objects.');

end

% Never do any enemy a small injury for they are like
% a snake which is half beaten; it will strike back 
% at the first chance it gets.
%
% Niccolo Machiavelli, "The Prince", Chapter 3

