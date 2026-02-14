% Refines a previously found step bracket by repeated interval
% bisection until a step satisfying Wolfe tests is found or the
% bracket collapses. Syntax:
%
%   [alpha,fx_1,gfx_1,exitflag,data]=...
%         sectioning(cost_function,A,B,x_0,fx_0,gfx_0,...
%                    dir,data,spin_system)
%
% Arguments:
%
%    cost_function     - objective function handle
%
%    a                 - lower bracket structure with
%                        fields alpha, fx, and gfx
%
%    b                 - upper bracket structure with
%                        fields alpha, fx, and gfx
%
%    x_0               - current optimisation vector
%
%    fx_0              - objective value at x_0
%
%    gfx_0             - gradient at x_0
%
%    dir               - search direction vector
%
%    data              - optimisation workspace structure
%
%    spin_system       - Spinach data structure with
%                        sectioning settings
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
% <https://spindynamics.org/wiki/index.php?title=sectioning.m>

function [alpha,fx_1,gfx_1,exitflag,data]=sectioning(cost_function,a,b,x_0,fx_0,...
                                                     gfx_0,dir,data,spin_system)

% Check consistency
grumble(cost_function,a,b,x_0,fx_0,gfx_0,dir);

% Apply coordinate freezing mask when requested
if ~isempty(spin_system.control.freeze)
    dir=dir.*(~spin_system.control.freeze(:));
    gfx_0=gfx_0.*(~spin_system.control.freeze(:));
end

% Iterate until Wolfe conditions 
% pass or the bracket collapses
while true

    % Simple interval bisection
    alpha=(a.alpha+b.alpha)/2;

    % Evaluate objective and gradient at current trial step
    [data,fx_1,gfx_1]=objeval(x_0+alpha*dir,cost_function,data,spin_system);

    % Apply coordinate freezing mask to the new gradient
    if ~isempty(spin_system.control.freeze)
        gfx_1=gfx_1.*(~spin_system.control.freeze(:));
    end

    % Save lower bracket endpoint
    lower_bracket_endpoint=a;

    % Update bracket based on Armijo and monotonicity tests
    if (~alpha_conds(1,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system))||...
       (~alpha_conds(0,[],a.fx,fx_1,[],[],[],[]))

        % Move upper bracket endpoint to current trial step
        b.alpha=alpha; b.fx=fx_1; b.gfx=gfx_1;

    else

        % Accept trial step if curvature condition is satisfied
        if alpha_conds(2,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system)
            exitflag=0; return;
        end

        % Move lower bracket endpoint to current trial step
        a.alpha=alpha; a.fx=fx_1; a.gfx=gfx_1;

        % Swap upper endpoint when derivative sign condition fails
        if (a.alpha-b.alpha)*(gfx_1'*dir)>=0
            b=lower_bracket_endpoint;
        end

    end

    % Terminate when bracket width becomes too small
    if abs((b.alpha-a.alpha)*(a.gfx'*dir))<eps(max(1,abs(a.fx)))

        % Current point is final
        alpha=a.alpha; fx_1=a.fx; gfx_1=a.gfx;

        % Keep improvements, otherwise terminate
        if alpha_conds(0,[],fx_0,fx_1,[],[],[],[])
            exitflag=0; return;
        else
            exitflag=-2; return;
        end

    end

end

end

% Consistency enforcement
function grumble(cost_function,a,b,x_0,fx_0,gfx_0,dir)
if ~isa(cost_function,'function_handle')
    error('cost_function must be a function handle.');
end
if ~isstruct(a)
    error('A must be a structure.');
end
if ~isstruct(b)
    error('B must be a structure.');
end
if (~isfield(a,'alpha'))||(~isfield(a,'fx'))||(~isfield(a,'gfx'))
    error('A must contain alpha, fx, and gfx fields.');
end
if (~isfield(b,'alpha'))||(~isfield(b,'fx'))||(~isfield(b,'gfx'))
    error('B must contain alpha, fx, and gfx fields.');
end
if isempty(a.alpha)||(~isnumeric(a.alpha))||(~isreal(a.alpha))||(~isscalar(a.alpha))
    error('A.alpha must be a real scalar.');
end
if isempty(b.alpha)||(~isnumeric(b.alpha))||(~isreal(b.alpha))||(~isscalar(b.alpha))
    error('B.alpha must be a real scalar.');
end
if isempty(a.fx)||(~isnumeric(a.fx))||(~isreal(a.fx))||(~isscalar(a.fx))
    error('A.fx must be a real scalar.');
end
if isempty(b.fx)||(~isnumeric(b.fx))||(~isreal(b.fx))||(~isscalar(b.fx))
    error('B.fx must be a real scalar.');
end
if isempty(a.gfx)||(~isnumeric(a.gfx))||(~isreal(a.gfx))||(~iscolumn(a.gfx))
    error('A.gfx must be a real column vector.');
end
if isempty(b.gfx)||(~isnumeric(b.gfx))||(~isreal(b.gfx))||(~iscolumn(b.gfx))
    error('B.gfx must be a real column vector.');
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
if isempty(dir)||(~isnumeric(dir))||(~isreal(dir))||(~iscolumn(dir))
    error('dir must be a real column vector.');
end
if ~isequal(size(dir),size(x_0),size(gfx_0),size(a.gfx),size(b.gfx))
    error('all vectors must have matching dimensions.');
end
end

% A gentleman a few rows in front of us took grave
% exception to the behaviour of an opposing player
% and identified him, very loudly, as the author 
% of the Critique of Pure Reason - repeatedly and
% with venom.
%
% Rod Liddle, in The Spectator

