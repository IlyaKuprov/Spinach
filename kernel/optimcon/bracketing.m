% Expands a trial step into a bracket that contains an acceptable
% line search point, or accepts the step directly if the Wolfe
% tests are met before sectioning becomes necessary.
%
% Syntax:
%
%    [A,B,alpha,fx,gfx,next_act,data]=...
%                bracketing(cost_function,alpha,dir,x_0,fx_0,...
%                           gfx_0,data,spin_system)
%
% Arguments:
%
%    cost_function     - objective function handle
%
%    alpha             - initial trial step length
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
%    spin_system       - Spinach data structure with
%                        line search settings
%
% Returns:
%
%    A                 - lower bracket point structure
%
%    B                 - upper bracket point structure
%
%    alpha             - accepted step length when found
%
%    fx                - objective value at accepted step
%
%    gfx               - gradient at accepted step
%
%    next_act          - continuation tag, either
%                        'sectioning' or 'none'
%
%    data              - updated optimisation workspace
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bracketing.m>

function [a,b,alpha,fx,gfx,next_act,data]=bracketing(cost_function,alpha,dir,x_0,...
                                                     fx_0,gfx_0,data,spin_system)

% Check consistency
grumble(cost_function,alpha,dir,x_0,fx_0,gfx_0);

% Apply coordinate freezing mask when requested
if ~isempty(spin_system.control.freeze)
    dir=dir.*(~spin_system.control.freeze(:));
    gfx_0=gfx_0.*(~spin_system.control.freeze(:));
end

% Evaluate objective and gradient at the first trial point
[data,fx_2,gfx_2]=objeval(x_0+alpha*dir,cost_function,data,spin_system);

% Apply coordinate freezing mask when requested
if ~isempty(spin_system.control.freeze)
    gfx_2=gfx_2.*(~spin_system.control.freeze(:));
end

% Initialise empty bracket records
a.alpha=[]; a.fx=[]; a.gfx=[];
b.alpha=[]; b.fx=[]; b.gfx=[];

% Initialise bracketing history variables
fx=fx_0; fx_1=fx_0;
gfx=gfx_0; gfx_1=gfx_0;
alpha_1=0; alpha_2=alpha;

% Expand bracket until acceptance or interval capture
while true

    % Capture bracket when Armijo or monotonicity fails
    if (~alpha_conds(1,alpha_2,fx_0,fx_2,gfx_0,[],dir,spin_system))||...
       (~alpha_conds(0,[],fx_1,fx_2,[],[],[],spin_system))

        % Store current interval endpoints
        a.alpha=alpha_1; a.fx=fx_1; a.gfx=gfx_1;
        b.alpha=alpha_2; b.fx=fx_2; b.gfx=gfx_2;

        % Hand over to sectioning stage
        next_act='sectioning'; return;

    end

    % Accept step immediately if curvature condition passes
    if alpha_conds(2,[],[],[],gfx_0,gfx_2,dir,spin_system)

        % Return accepted step information
        alpha=alpha_2; fx=fx_2; gfx=gfx_2;

        % No further line search steps are required
        next_act='none'; return;

    end

    % Capture bracket when directional derivative changes sign
    if ~alpha_conds(3,[],[],[],[],gfx_2,dir,spin_system)

        % Store current interval endpoints
        a.alpha=alpha_2; a.fx=fx_2; a.gfx=gfx_2;
        b.alpha=alpha_1; b.fx=fx_1; b.gfx=gfx_1;

        % Hand over to sectioning stage
        next_act='sectioning'; return;

    end

    % Build interpolation window ahead of current trial point
    br_end_pt_A=2*alpha_2-alpha_1;
    br_end_pt_B=alpha_2+spin_system.control.ls_tau1*(alpha_2-alpha_1);

    % Maximise the cubic model inside interpolation bounds
    alpha_new=cubic_interp(br_end_pt_A,br_end_pt_B,alpha_1,alpha_2,...
                           fx_1,gfx_1'*dir,fx_2,gfx_2'*dir);

    % Shift history to the new trial point
    alpha_1=alpha_2; alpha_2=alpha_new;
    fx_1=fx_2; gfx_1=gfx_2; x_1=x_0+alpha_2*dir;

    % Evaluate objective and gradient at the new trial point
    [data,fx_2,gfx_2]=objeval(x_1,cost_function,data,spin_system);

    % Apply coordinate freezing mask to the new gradient
    if ~isempty(spin_system.control.freeze)
        gfx_2=gfx_2.*(~spin_system.control.freeze(:));
    end

end

end

% Consistency enforcement
function grumble(cost_function,alpha,dir,x_0,fx_0,gfx_0)
if ~isa(cost_function,'function_handle')
    error('cost_function must be a function handle.');
end
if isempty(alpha)||(~isnumeric(alpha))||(~isreal(alpha))||(~isscalar(alpha))
    error('alpha must be a real scalar.');
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

% Одной рукой бунтую, другой пишу донос.
%
% Михаил Щербаков

