% L1 norm Tikhonov regularised solver for A*x=y where
% A is an ill-conditioned matrix. The error functional
% is norm(A*x-y,2)^2+lambda*norm(x,1), it is minimised
% using the FISTA algorithm. The user specifies the de-
% sired number of non-zeroes, lambda parameter is then
% found by bracketing / bisection. Syntax:
%
%           [x,err,reg]=tikhol1n(A,y,nnzt)
%
% Parameters:
%
%    A    - a real or complex matrix
%
%    y    - a real or complex column vector
%
%    nnzt - the target for the number of 
%           non-zeroes in the solution
%
% Outputs:
%
%    x   - a real or complex vector
%
%    err - squared 2-norm of the fitting
%          error divided by the squared
%          2-norm of the solution
%
%    reg - 1-norm of the solution
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=tikhol1n.m>

function [x,err,reg]=tikhol1n(A,y,nnzt)

% Check consistency
grumble(A,y,nnzt);

% Tolerances
normest_tol=1e-3;   % relative 2-norm estimation tolerance
step_norm_tol=1e-6; % relative step norm convergence tolerance

% Pre-compute CT
A_ct=ctranspose(A);
        
% Lipschitz constant and initial threshold
L=2*(1+normest_tol)*normest(A,normest_tol)^2; thr=1/L;

% Complex soft thresholding function
soft_thr=@(x,thr)sign(x).*max(abs(x)-thr,0);
        
% Zero initial guess
x=zeros(size(A,2),1); x_old=x;

% Threshold brackets for ZF targeting
thr_lower=0; thr_upper=max(abs(A_ct*y));

% Iteration counters and momentum
t=1; iter_count=0; converged=false;

% Tell the user we are starting to iterate
disp(['FISTA called with target nnz = ' int2str(nnzt)]);
        
% FISTA loop
while ~converged

    % Error vector
    err_vec=A*x-y;

    % Error gradient 
    g=2*(A_ct*err_vec);

    % Compute proximal point
    x_prox=soft_thr(x-g/L,thr);

    % Update Nesterov momentum
    t_new=0.5*(1+sqrt(1+4*t^2));

    % Take the step
    x=x_prox+((t-1)/t_new)*(x_prox-x_old);

    % Check convergence
    soln_norm=norm(x_prox,2);
    step_norm=norm(x_prox-x_old,2);
    converged=(step_norm<step_norm_tol*soln_norm);
    
    % Close the loop
    t=t_new; x_old=x_prox;
    iter_count=iter_count+1;

    % Progress report and nnz targeting
    if mod(iter_count,1000)==0 || converged

        % Check solution health and the nnz target within ±1 tolerance
        if nnz(x)==0 || (converged && ~ismember(nnz(x),[nnzt-1, nnzt, nnzt+1]))

            % Decisions
            if nnz(x)>nnzt

                % Need fewer zeroes  
                thr_lower=thr;

            else
                
                % Need more zeroes
                thr_upper=thr;

            end
            
            % Recalculate midpoint
            thr=(thr_lower+thr_upper)/2;

            % Restart the calculation
            t=1; converged=false;
            
            % Inform the user
            disp(['zero threshold brackets updated: lower ' ...
                  num2str(thr_lower) ', upper ' num2str(thr_upper)]);
            
            % Detect stagnation
            if (thr_upper-thr_lower)/(thr_upper+thr_lower)<1e-6
                error('nnz target unreachable');
            end

        else

            % Get solver state metrics
            err=norm(err_vec,2)^2/norm(x,2)^2; reg=norm(x,1);

            % Print the report
            disp(['FISTA iter ' int2str(iter_count) ', nnz ' int2str(nnz(x)) ...
                  ', rel. sq. err. ' num2str(err) ', 1-norm ' num2str(reg) ...
                  ', rel. step ' num2str(step_norm/soln_norm)                ...
                  ', zero thr. ' num2str(thr)]);

        end

    end

end

end

% Consistency enforcement
function grumble(A,y,nnzt)
if (~isnumeric(A))||(~ismatrix(A))
    error('A must be a numeric matrix.');
end
if (~isnumeric(y))||(~iscolumn(y))
    error('y must be a numeric column vector.');
end
if size(A,1)~=size(y,1)
    error('the number of rows in A must match the number of elements in y.');
end
if (~isnumeric(nnzt))||(~isreal(nnzt))||(~isscalar(nnzt))||...
   (nnzt<1)||(mod(nnzt,1)~=0)
    error('nnzt must be a positive integer.');
end
if nnzt>size(A,2)
    error('nnzt cannot exceed the number of columns in A.');
end
end

% "He begins working on calculus problems in his head as 
% soon as he awakens. He did calculus while driving in his
% car, while sitting in the living room, and while lying 
% in bed at night."
%
% Mary Louise Bell (Richard Feynman's second
% wife), in her divorce complaint

