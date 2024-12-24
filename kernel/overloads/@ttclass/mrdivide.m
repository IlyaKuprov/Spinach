% Divides a tensor train object by a scalar. Syntax:
%
%                    c=mrdivide(a,b)
%
% The first input should be a ttclass object and the 
% second one a scalar.
%
% d.savostyanov@soton.ac.uk
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/mrdivide.m>

function a=mrdivide(a,b)

% Division of tensor train by a scalar
if isa(a,'ttclass')&&isscalar(b)

    % Divide the coefficients and update the tolerances
    a.coeff=a.coeff/b; a.tolerance=a.tolerance/abs(b);
    
else
    
    % Complain and bomb out
    error('the first argument must be a tensor train and the second one a scalar.');
    
end
    
end

% It is dangerous to be right in matters on which the established
% authorities are wrong. 
%
% Voltaire

