% Refines a previously found step bracket by repeated cubic
% interpolation until a step satisfying Wolfe tests is found
% or the bracket collapses below numerical resolution.
%
% Syntax:
%
%   [alpha,fx_1,gfx_1,exitflag,data]=sectioning(cost_function,...
%                      A,B,x_0,fx_0,gfx_0,dir,data,spin_system)
%
% Arguments:
%
%    cost_function     - objective function handle
%
%    A                 - lower bracket structure with
%                        fields alpha, fx, and gfx
%
%    B                 - upper bracket structure with
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

function [alpha,fx_1,gfx_1,exitflag,data]=sectioning(cost_function,A,B,x_0,fx_0,gfx_0,dir,data,spin_system)

% Check consistency
grumble(cost_function,A,B,x_0,fx_0,gfx_0,dir,data,spin_system);

% Apply coordinate freezing mask when requested
if ~isempty(spin_system.control.freeze)
    dir=dir.*(~spin_system.control.freeze(:));
    gfx_0=gfx_0.*(~spin_system.control.freeze(:));
end

% Iterate until Wolfe conditions pass or bracket collapses
while true

    % Build reduced interpolation bounds inside the bracket
    end_A=A.alpha+spin_system.control.ls_tau2*(B.alpha-A.alpha);
    end_B=B.alpha-spin_system.control.ls_tau3*(B.alpha-A.alpha);

    % Maximise cubic model inside reduced interpolation bounds
    alpha=cubic_interp(end_A,end_B,A.alpha,B.alpha,A.fx,A.gfx'*dir,B.fx,B.gfx'*dir);

    % Stop when interpolation displacement is numerically unresolved
    if abs((alpha-A.alpha)*(A.gfx'*dir))<=eps(max(1,abs(fx_0)))
        exitflag=-2; return;
    end

    % Evaluate objective and gradient at current trial step
    [data,fx_1,gfx_1]=objeval(x_0+alpha*dir,cost_function,data,spin_system);

    % Apply coordinate freezing mask to the new gradient
    if ~isempty(spin_system.control.freeze)
        gfx_1=gfx_1.*(~spin_system.control.freeze(:));
    end

    % Store current lower bracket endpoint before updates
    tmp=A;

    % Update bracket based on Armijo and monotonicity tests
    if (~alpha_conds(1,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system))||...
       (~alpha_conds(0,alpha,A.fx,fx_1,A.gfx,gfx_1,dir,spin_system))

        % Move upper bracket endpoint to current trial step
        B.alpha=alpha; B.fx=fx_1; B.gfx=gfx_1;

    else

        % Accept trial step if curvature condition is satisfied
        if alpha_conds(2,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system)
            exitflag=0; return;
        end

        % Move lower bracket endpoint to current trial step
        A.alpha=alpha; A.fx=fx_1; A.gfx=gfx_1;

        % Swap upper endpoint when derivative sign condition fails
        if (A.alpha-B.alpha)*(gfx_1'*dir)>=0
            B=tmp;
        end

    end

    % Terminate when bracket width falls below machine precision
    if abs(B.alpha-A.alpha)<eps
        alpha=A.alpha; fx_1=A.fx;
        gfx_1=A.gfx; exitflag=-2; return;
    end

end

end

% Consistency enforcement
function grumble(cost_function,A,B,x_0,fx_0,gfx_0,dir,data,spin_system)
if ~isa(cost_function,'function_handle')
    error('cost_function must be a function handle.');
end
if ~isstruct(A)
    error('A must be a structure.');
end
if ~isstruct(B)
    error('B must be a structure.');
end
if (~isfield(A,'alpha'))||(~isfield(A,'fx'))||(~isfield(A,'gfx'))
    error('A must contain alpha, fx, and gfx fields.');
end
if (~isfield(B,'alpha'))||(~isfield(B,'fx'))||(~isfield(B,'gfx'))
    error('B must contain alpha, fx, and gfx fields.');
end
if isempty(A.alpha)||(~isnumeric(A.alpha))||(~isreal(A.alpha))||(~isscalar(A.alpha))
    error('A.alpha must be a real scalar.');
end
if isempty(B.alpha)||(~isnumeric(B.alpha))||(~isreal(B.alpha))||(~isscalar(B.alpha))
    error('B.alpha must be a real scalar.');
end
if isempty(A.fx)||(~isnumeric(A.fx))||(~isreal(A.fx))||(~isscalar(A.fx))
    error('A.fx must be a real scalar.');
end
if isempty(B.fx)||(~isnumeric(B.fx))||(~isreal(B.fx))||(~isscalar(B.fx))
    error('B.fx must be a real scalar.');
end
if isempty(A.gfx)||(~isnumeric(A.gfx))||(~isreal(A.gfx))||(~iscolumn(A.gfx))
    error('A.gfx must be a real column vector.');
end
if isempty(B.gfx)||(~isnumeric(B.gfx))||(~isreal(B.gfx))||(~iscolumn(B.gfx))
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
if ~isequal(size(dir),size(x_0),size(gfx_0),size(A.gfx),size(B.gfx))
    error('all vectors must have matching dimensions.');
end
if ~isstruct(data)
    error('data must be a structure.');
end
if (~isstruct(spin_system))||(~isfield(spin_system,'control'))
    error('spin_system.control must be present.');
end
if (~isfield(spin_system.control,'ls_tau2'))||(~isnumeric(spin_system.control.ls_tau2))||...
   (~isreal(spin_system.control.ls_tau2))||(~isscalar(spin_system.control.ls_tau2))
    error('spin_system.control.ls_tau2 must be a real scalar.');
end
if (~isfield(spin_system.control,'ls_tau3'))||(~isnumeric(spin_system.control.ls_tau3))||...
   (~isreal(spin_system.control.ls_tau3))||(~isscalar(spin_system.control.ls_tau3))
    error('spin_system.control.ls_tau3 must be a real scalar.');
end
if ~isfield(spin_system.control,'freeze')
    error('spin_system.control.freeze must be present.');
end
if ~isempty(spin_system.control.freeze)
    if (~islogical(spin_system.control.freeze))||...
       (~isvector(spin_system.control.freeze))||...
       (numel(spin_system.control.freeze)~=numel(dir))
        error('spin_system.control.freeze must be a logical vector matching dir length.');
    end
end
end


