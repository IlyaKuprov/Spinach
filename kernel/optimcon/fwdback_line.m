% Locates an acceptable ascent step using a simple forward- and
% backtracking line search strategy. Syntax:
%
%    [alpha,fx_1,gfx_1,exitflag,data]=...
%        fwdback_line(cost_function,alpha_0,dir,x_0,fx_0,...
%                     gfx_0,data,spin_system)
%
% Arguments:
%
%    cost_function     - objective function handle
%
%    alpha_0           - initial trial step length
%
%    dir               - search direction vector
%
%    x_0               - current optimisation vector
%
%    fx_0              - objective value at x_0
%
%    gfx_0             - gradient at x_0
%
%    data              - optimisation workspace structure
%
%    spin_system       - Spinach data structure with line
%                        search settings in control
%
% Returns:
%
%    alpha             - accepted step length
%
%    fx_1              - objective value at alpha
%
%    gfx_1             - gradient at alpha
%
%    exitflag          - 0 on success, -2 on failure
%
%    data              - updated optimisation workspace
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=fwdback_line.m>

function [alpha,fx_1,gfx_1,exitflag,data]=fwdback_line(cost_function,alpha_0,dir,x_0,...
                                                       fx_0,gfx_0,data,spin_system)

% Check consistency
grumble(cost_function,alpha_0,dir,x_0,fx_0,gfx_0);

% Set line search defaults
exitflag=-2;
alpha=alpha_0;
fx_1=fx_0;
gfx_1=gfx_0;
max_trials=25;

% Apply coordinate freezing mask when requested
if ~isempty(spin_system.control.freeze)
    dir=dir.*(~spin_system.control.freeze(:));
    gfx_0=gfx_0.*(~spin_system.control.freeze(:));
end

% Evaluate the first trial step
[data,fx_trial,gfx_trial]=objeval(x_0+alpha*dir,cost_function,data,spin_system);

% Apply coordinate freezing mask to the trial gradient
if ~isempty(spin_system.control.freeze)
    gfx_trial=gfx_trial.*(~spin_system.control.freeze(:));
end

% Decide between forward- and backtracking branches
if alpha_conds(1,alpha,fx_0,fx_trial,gfx_0,gfx_trial,dir,spin_system)&&...
   alpha_conds(0,[],fx_0,fx_trial,[],[],[],spin_system)

    % Accept the first trial point as current best
    exitflag=0;
    fx_1=fx_trial;
    gfx_1=gfx_trial;

    % Expand trial step while objective continues improving
    for n=1:max_trials

        % Propose a larger trial step
        alpha_try=spin_system.control.ls_tau1*alpha;

        % Evaluate objective and gradient at the expanded point
        [data,fx_try,gfx_try]=objeval(x_0+alpha_try*dir,cost_function,data,spin_system);

        % Apply coordinate freezing mask to the expanded gradient
        if ~isempty(spin_system.control.freeze)
            gfx_try=gfx_try.*(~spin_system.control.freeze(:));
        end

        % Stop forward tracking when sufficient increase fails
        if ~alpha_conds(1,alpha_try,fx_0,fx_try,gfx_0,gfx_try,dir,spin_system)
            break;
        end

        % Stop forward tracking when monotonic growth fails
        if ~alpha_conds(0,[],fx_1,fx_try,[],[],[],spin_system)
            break;
        end

        % Accept the expanded point and continue exploring
        alpha=alpha_try;
        fx_1=fx_try;
        gfx_1=gfx_try;

    end

else

    % Contract the step length until sufficient increase is reached
    for n=1:max_trials

        % Contract trial step length
        alpha=spin_system.control.ls_tau3*alpha;

        % Stop when trial step falls below numerical resolution
        if alpha<eps
            return;
        end

        % Evaluate objective and gradient at the contracted point
        [data,fx_trial,gfx_trial]=objeval(x_0+alpha*dir,cost_function,data,spin_system);

        % Apply coordinate freezing mask to the contracted gradient
        if ~isempty(spin_system.control.freeze)
            gfx_trial=gfx_trial.*(~spin_system.control.freeze(:));
        end

        % Accept point when sufficient increase and monotonicity pass
        if alpha_conds(1,alpha,fx_0,fx_trial,gfx_0,gfx_trial,dir,spin_system)&&...
           alpha_conds(0,[],fx_0,fx_trial,[],[],[],spin_system)

            % Store accepted trial point
            exitflag=0;
            fx_1=fx_trial;
            gfx_1=gfx_trial;
            return;

        end

    end

end

end

% Consistency enforcement
function grumble(cost_function,alpha_0,dir,x_0,fx_0,gfx_0)
if ~isa(cost_function,'function_handle')
    error('cost_function must be a function handle.');
end
if isempty(alpha_0)||(~isnumeric(alpha_0))||(~isreal(alpha_0))||(~isscalar(alpha_0))
    error('alpha_0 must be a real scalar.');
end
if isempty(dir)||(~isnumeric(dir))||(~isreal(dir))||(~iscolumn(dir))
    error('dir must be a real column vector.');
end
if isempty(x_0)||(~isnumeric(x_0))||(~isreal(x_0))||(~iscolumn(x_0))
    error('x_0 must be a real column vector.');
end
if isempty(fx_0)||(~isnumeric(fx_0))||(~isreal(fx_0))||(~isscalar(fx_0))
    error('fx_0 must be a real scalar.');
end
if isempty(gfx_0)||(~isnumeric(gfx_0))||(~isreal(gfx_0))||(~iscolumn(gfx_0))
    error('gfx_0 must be a real column vector.');
end
if ~isequal(size(dir),size(x_0),size(gfx_0))
    error('dir, x_0, and gfx_0 must have matching dimensions.');
end
end


