% Estimate of the Jacobian matrix of a vector valued 
% function of n variables. Syntax:
%
%              [jac,err] = jacobianest(fun,x0)
%
% Parameters:
%
%  fun - (vector valued) analytical function to differentiate.
%        fun must be a function of the vector or array x0.
% 
%  x0  - vector location at which to differentiate fun
%        If x0 is an nxm array, then fun is assumed to be
%        a function of n*m variables.
%
% Outputs:
%
%  jac - array of first partial derivatives of fun.
%        Assuming that x0 is a vector of length p
%        and fun returns a vector of length n, then
%        jac will be an array of size (n,p)
%
%  err - vector of error estimates corresponding to
%        each partial derivative in jac.
%
% John D'Errico
%
% <https://spindynamics.org/wiki/index.php?title=jacobianest.m>

function [jac,err] = jacobianest(fun,x0)

% Check consistency
grumble(fun,x0);

% get the length of x0 for the size of jac
nx = numel(x0);

% Finite difference gridding
MaxStep = 100;
StepRatio = 2.0000001;
relativedelta = MaxStep*StepRatio .^(0:-1:-25);
nsteps = length(relativedelta);

% was a string supplied?
if ischar(fun)
  fun = str2func(fun);
end

% get fun at the center point
f0 = fun(x0);
f0 = f0(:);
n = length(f0);
if n==0
  % empty begets empty
  jac = zeros(0,nx);
  err = jac;
  return
end

% total number of derivatives we will need to take
jac = zeros(n,nx);
err = jac;
for i = 1:nx
  x0_i = x0(i);
  if x0_i ~= 0
    delta = x0_i*relativedelta;
  else
    delta = relativedelta;
  end
  
  % evaluate at each step, centered around x0_i
  % difference to give a second order estimate
  fdel = zeros(n,nsteps);
  for j = 1:nsteps
    fdif = fun(swapelement(x0,i,x0_i + delta(j))) - fun(swapelement(x0,i,x0_i - delta(j)));
    
    fdel(:,j) = fdif(:);
  end
  
  % these are pure second order estimates of the
  % first derivative, for each trial delta.
  derest = fdel.*repmat(0.5 ./ delta,n,1);
  
  % The error term on these estimates has a second order
  % component, but also some 4th and 6th order terms in it.
  % Use Romberg exrapolation to improve the estimates to
  % 6th order, as well as to provide the error estimate.
  
  % loop here, as rombextrap coupled with the trimming
  % will get complicated otherwise.
  for j = 1:n
    [der_romb,errest] = rombextrap(StepRatio,derest(j,:),[2 4]);
    
    % trim off 3 estimates at each end of the scale
    nest = length(der_romb);
    trim = [1:3, nest+(-2:0)];
    [der_romb,tags] = sort(der_romb);
    der_romb(trim) = [];
    tags(trim) = [];
    errest = errest(tags);
    
    % now pick the estimate with the lowest predicted error
    [err(j,i),ind] = min(errest);
    jac(j,i) = der_romb(ind);
  end
  
end

end

% Swaps val as element ind, into the vector vec
function vec = swapelement(vec,ind,val)
    vec(ind) = val;
end 

% Romberg extrapolation
%
% StepRatio - Ratio decrease in step
% der_init - initial derivative estimates
% rombexpon - higher order terms to cancel using the romberg step
%
% der_romb - derivative estimates returned
% errest - error estimates
% amp - noise amplification factor due to the romberg step

function [der_romb,errest] = rombextrap(StepRatio,der_init,rombexpon)

% Inverse step ratio
srinv = 1/StepRatio;

% do nothing if no romberg terms
nexpon = length(rombexpon);
rmat = ones(nexpon+2,nexpon+1);

% two romberg terms
rmat(2,2:3) = srinv.^rombexpon;
rmat(3,2:3) = srinv.^(2*rombexpon);
rmat(4,2:3) = srinv.^(3*rombexpon);

% qr factorization used for the extrapolation as well
% as the uncertainty estimates
[qromb,rromb] = qr(rmat,0);

% this does the extrapolation to a zero step size.
ne = length(der_init);
rhs = vec2mat(der_init,nexpon+2,ne - (nexpon+2));
rombcoefs = rromb\(qromb'*rhs);
der_romb = rombcoefs(1,:)';

% uncertainty estimate of derivative prediction
s = sqrt(sum((rhs - rmat*rombcoefs).^2,1));
rinv = rromb\eye(nexpon+1);
cov1 = sum(rinv.^2,2); % 1 spare dof
errest = s'*12.7062047361747*sqrt(cov1(1));

end

% Forms the matrix M, such that M(i,j) = vec(i+j-1)
function mat = vec2mat(vec,n,m)
[i,j] = ndgrid(1:n,0:m-1);
ind = i+j; mat = vec(ind);
if n==1, mat = mat'; end
end

% Consistency enforcement
function grumble(fun,x0)
if ~isa(fun,'function_handle')
    error('fun must be a function handle.');
end
if ~isnumeric(x0)
    error('x0 must be numeric.');
end
end

% To prove that you are not a robot, injure a human 
% being, or, through inaction, allow a human being
% to come to harm.

