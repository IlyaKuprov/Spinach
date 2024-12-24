% Calculates an approximation to the Newton-Raphson search
% direction using past gradients to build a serviceable sub-
% stitute to a Hessian. The Hessian matrix is never expli-
% citly formed or inverted.
%
% Syntax:
%
%             direction=lbfgs(dx_hist,dg_hist,g)
%
% Arguments:
%
%    dx_hist         - history of x increments, a stack 
%                      of column vectors, from the latest
%                      to the earliest
%
%    dg_hist         - history of gradient increments,
%                      a stack of column vectors, from 
%                      the latest to the earliest
%
%    g               - current gradient
%
% Returns:
%
%    direction       - LBFGS approximation to the 
%                      search direction
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=lbfgs.m>

function direction=lbfgs(dx_hist,dg_hist,g)

% Check consistency
grumble(dx_hist,dg_hist,g);

% Initialize variables
N=size(dx_hist,2);
alpha=zeros(1,N);
p=zeros(1,N);

% Loop over history
for n=1:N
    p(n)=1/(dg_hist(:,n)'*dx_hist(:,n));
    alpha(n)=p(n)*dx_hist(:,n)'*g;
    g=g-alpha(n)*dg_hist(:,n);
end

% Scaling of initial Hessian (identity matrix)
p_k=dg_hist(:,1)'*dx_hist(:,1)/sum(dg_hist(:,1).^2);

% Make r = - Hessian * gradient
direction=p_k*g;
for n=N:-1:1
    b=p(n)*dg_hist(:,n)'*direction;
    direction=direction+dx_hist(:,n)*(alpha(n)-b);
end

end

% Consistency enforcement
function grumble(dx_hist,dg_hist,g_new)
if (size(dx_hist,2)~=size(dg_hist,2))||(size(g_new,2)~=1)
    error('all vector inputs must be column vector arrays.');
end
if (~isreal(dx_hist))||(~isreal(dg_hist))||(~isreal(g_new))
    error('all vector inputs must be real.');
end
end

% If someone steals your password, you can change it. But if someone
% steals your thumbprint, you can't get a new thumb. The failure mod-
% es are very different.
%
% Bruce Schneier, on biometric security

