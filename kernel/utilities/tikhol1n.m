% L1 norm Tikhonov regularised solver for A*x=y where
% A is an ill-conditioned matrix. The error functional
% is norm(A*x-y,1)^2+lambda*norm(x,1), it is minimised
% using the FISTA algorithm. Syntax:
%
%           [x,err,reg]=tikhol1n(A,y,zft)
%
% Parameters:
%
%    A - a real or complex matrix
%
%    y - a real or complex column vector
%
%    zft - zero fraction target
%
% Outputs:
%
%    x - a real or complex vector
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=tikhol1n.m>

function [x,err,reg]=tikhol1n(A,y,zft)

% Tolerances
normest_tol=1e-3;   % relative 2-norm estimation tolerance
step_norm_tol=1e-6; % relative step norm convergence tolerance

% Pre-compute CT
A_ct=ctranspose(A);
        
% Lipschitz constant and initial threshold
L=2*normest(A,normest_tol)^2; thr=1/L;

% Complex soft thresholding function
soft_thr=@(x,thr)sign(x).*max(abs(x)-thr,0);
        
% Zero initial guess
x=zeros(size(A,2),1); x_old=x;

% Iteration counters and momentum
t=1; iter_count=0; converged=false;

% Report to the user
disp(['FISTA called with zero fraction target = ' num2str(zft)]);
        
% FISTA iteration loop
while ~converged

    % Get the error vector
    err_vec=A*x-y;

    % Get the gradient of the error
    g=2*(A_ct*err_vec);

    % Compute the proximal point
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

    % Analyse and report progress
    if mod(iter_count,1000)==0 || converged

        % Get the zero fraction
        zf=1-nnz(x)/numel(x);

        % Check the zero fraction target
        if abs(zft-zf)>0.01 && converged

            % Update the threshold
            thr=thr*2^sign(zft-zf);

            % Reset the state
            t=1; converged=false;

            % Let the user know
            disp('threshold updated');

        else

            % Get solver state metrics
            err=norm(err_vec,2)^2; reg=norm(x,1);

            % Print the report
            disp(['FISTA iter ' int2str(iter_count) ', zf ' num2str(zf) ...
                  ', err '  num2str(err) ', reg ' num2str(reg)          ...
                  ', step ' num2str(step_norm/soln_norm)                ...
                  ', thr ' num2str(thr)]);

        end

    end

end

% Compute figures of merit
err=norm(err_vec,2)^2; reg=norm(x,1);

end

% "He begins working on calculus problems in his head as 
% soon as he awakens. He did calculus while driving in his
% car, while sitting in the living room, and while lying 
% in bed at night."
%
% Mary Louise Bell (Richard Feynman's second
% wife), in her divorce complaint

