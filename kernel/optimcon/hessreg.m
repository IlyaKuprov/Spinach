% RFO regularisation for Newton-Raphson Hessian and gradient 
% pairs. Syntax:
%
%           [H,data]=hessreg(spin_system,H,g,data)
%
% Parameters:
%
%   H    - Hessian matrix to be regularised
%
%   g    - gradient computed at the same point as H
%
%   data - diagnostic data structure
%
% Outputs:
%
%   H    - regularised Hessian
%
%   data - updated diagnostic data structure with
%          data.count.rfo incremented by the number
%          of RFO iterations taken
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=hessreg.m>

function [H,data]=hessreg(spin_system,H,g,data)

% Check consistency
grumble(H,g);

% Set shorthands
alpha=spin_system.control.reg_alpha;
phi=spin_system.control.reg_phi;
max_iter=spin_system.control.reg_max_iter;
max_cond=spin_system.control.reg_max_cond;

% Check if RFO is needed
[~,p]=chol(H); pos_def=~logical(p);
if pos_def&&(cond(H,2)<max_cond)
    return;
end

% Start RFO iteration loop
for ind=1:max_iter
    
    % Make auxiliary Hessian
    H=[(alpha^2)*H   alpha*g;
        alpha*g'     0];
    
    % Calculate eigenvalue shift
    sigma=min([0 min(eig(H,'vector'))]);
   
    % Apply eigenvalue shift
    H=H-sigma*speye(size(H));
    
    % Truncate and scale back
    H=H(1:(end-1),(1:end-1));
    H=H/(alpha^2);
    
    % Update alpha
    alpha=alpha*phi;
    
    % Increment the counter
    data.count.rfo=data.count.rfo+1;
    
    % Break if condition number reached
    if cond(H,2)<max_cond, break; end
    
end

% Clean up the result
H=real(H+H')/2;

end

% Consistency enforcement
function grumble(H,g)
if (~isnumeric(H))||(size(H,1)~=size(H,2))||...
   (~issymmetric(H))||(~isreal(H))
    error('H must be a real square symmatric matrix.');
end
if (~isnumeric(g))||(~isreal(g))||(~iscolumn(g))
    error('g must be a real column vector.');
end
if (numel(g)~=size(H,1))
    error('dimensions of g and H must match.');
end
end

% Единственное, что я понимаю в арбузах - это если я по 
% арбузу постучала, а из него постучали в ответ, то вот
% это прям однозначно плохой арбуз.
%
% Russian internet folklore

