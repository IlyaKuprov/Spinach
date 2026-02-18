% Calculates a BFGS approximation to the Newton-Raphson search
% direction for maximising a function using past gradients to
% build a serviceable substitute to a Hessian. Unlike LBFGS, 
% the pseudo-Hessian matrix is formed explicitly. Syntax:
%
%                 H=bfgs(dx_hist,dg_hist,g)
%
% Arguments:
%
%    dx_hist  - history of x increments, a stack
%               of column vectors, from the latest
%               to the earliest
%
%    dg_hist  - history of gradient increments,
%               a stack of column vectors, from
%               the latest to the earliest
%
%    g        - current gradient (used for sizing)
%
% Returns:
%
%    H        - BFGS approximation to the Hessian 
%               matrix corresponding to the *nega-
%               tive* Hessian of the objective.
%               The corresponding ascent directi-
%               on is obtained as:  direction=H\g
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bfgs.m>

function H=bfgs(dx_hist,dg_hist,g)

% Check consistency
grumble(dx_hist,dg_hist,g);

% Curvature tolerance
curv_rel_tol=0.01;

% Relevant inner products
dgdx=sum(dg_hist.*dx_hist,1);  
dgdg=sum(dg_hist.*dg_hist,1);  
dxdx=sum(dx_hist.*dx_hist,1);   

% Curvature pair validation
valid=isfinite(dgdx)&isfinite(dgdg)&isfinite(dxdx)&...
      (dgdg>0)&(dxdx>0)&(dgdx<-curv_rel_tol*sqrt(dgdg.*dxdx));

% Dropping bad pairs
dx_hist=dx_hist(:,valid);
dg_hist=dg_hist(:,valid);

% Return identity if 
% all pairs are bad
if isempty(dx_hist)
    H=eye(numel(g)); return;
end

% Dimension and history length
n=numel(g); N=size(dx_hist,2);

% Pull the first curvature pair
dx=dx_hist(:,1); dg=-dg_hist(:,1);

% Inner products
dgdx=dg'*dx; dgdg=dg'*dg;

% Unit matrix scaling
if (~isfinite(dgdx))||...
   (~isfinite(dgdg))||...
   (dgdx<=eps)||(dgdg<=0)

    % Bad pair
    gamma=1;

else

    % Good pair
    gamma=dgdg/dgdx;

end

% Initial guess
H=gamma*eye(n);

% BFGS updates
for k=N:-1:1

    % Pull next curvature pair
    dx=dx_hist(:,k); dg=-dg_hist(:,k);

    % Make BFGS components
    dgdx=dg'*dx; Hdx=H*dx; dxHdx=dx'*Hdx;

    % Run pair validity checks to protect numerics
    if (~isfinite(dgdx))||(dgdx<=eps), continue; end
    if (~isfinite(dxHdx))||(dxHdx<=eps), continue; end

    % Plain BFGS Hessian update
    H=H-(Hdx*Hdx')/dxHdx+(dg*dg')/dgdx;

    % Numerical clean-up
    H=real((H+H')/2);

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

% Of all tyrannies, a tyranny sincerely exercised for the good of its
% victims may be the most oppressive. It would be better to live under
% robber barons than under omnipotent moral busybodies. The robber ba-
% ron's cruelty may sometimes sleep, his cupidity may at some point be
% satiated; but those who torment us for our own good will torment us
% without end for they do so with the approval of their own conscience.
%
% C.S. Lewis

