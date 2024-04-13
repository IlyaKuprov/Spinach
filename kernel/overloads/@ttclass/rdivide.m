% Divides a tensor train object by a scalar. Syntax:
%
%                    c=rdivide(a,b)
%
% The first input should be a ttclass object and the 
% second one a scalar.
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/rdivide.m>

function a=rdivide(a,b)

% Division of tensor train by a scalar
if isa(a,'ttclass')&&isscalar(b)

    % Divide the coefficients and update the tolerances
    a.coeff=a.coeff/b; a.tolerance=a.tolerance/abs(b);
    
else
    
    % Complain and bomb out
    error('the first argument must be a tensor train and the second one a scalar.');
    
end
    
end

% Documentation is like sex: when it is good, it is 
% very, very good, and when it is bad it's still bet-
% ter than nothing.
%
% Jim Hargrove

