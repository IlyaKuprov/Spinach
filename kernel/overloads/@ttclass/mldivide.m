% Solves a linear system with tensor train objects. Syntax:
%
%                      x=mldivide(A,y)
%
% Parameters:
%
%    A - ttclass matrix 
%
%    x - ttclass vector
%
% Outputs:
%
%    x - ttclass vector
%
% Note: the AMEn-solve algorithm is applied to symmetrised 
%       system (A'*A)*x=A'*y
%
% d.savostyanov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/mldivide.m>

function x=mldivide(A,y)

if isa(A,'ttclass')&&isa(y,'ttclass')
    
    % Shrink the operands
    A=shrink(A); 
    y=shrink(y);
    
    % Form a symmetrised system
    AA=shrink(A'*A); 
    Ay=shrink(A'*y);
    
    % Solve it with AMEn algorithm
    x=amensolve(AA,Ay,1e-6);
    
else
    
    % Complain and bomb out
    error('both arguments should be tensor trains.');
    
end

end

% The penalty for success is to be bored by the people
% who used to snub you.
%
% Nancy Astor

