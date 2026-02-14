% Applies one of the line search acceptance tests used by the brac-
% keting and sectioning routines in constrained optimisation and 
% returns true when the chosen condition is satisfied. Syntax:
%
%          test=alpha_conds(test_type,alpha,fx_0,fx_1,...
%                           gfx_0,gfx_1,dir,spin_system)
%
% Arguments:
%
%    test_type         - condition selector:
%
%                        0 for monotonic increase test
%                        1 for Armijo sufficient increase test
%                        2 for strong Wolfe curvature test
%                        3 for ascent direction test
%
%    alpha             - trial step length
%
%    fx_0              - objective value at the initial point
%
%    fx_1              - objective value at the trial point
%
%    gfx_0             - gradient at the initial point
%
%    gfx_1             - gradient at the trial point
%
%    dir               - search direction vector
%
%    spin_system       - Spinach data structure with line
%                        search settings in control
%
% Returns:
%
%    test              - logical true if the selected
%                        condition is satisfied
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=alpha_conds.m>

function test=alpha_conds(test_type,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system)

% Check consistency
grumble(test_type,alpha,fx_0,fx_1,gfx_0,gfx_1,dir);

% Test selection
if test_type==0

    % Monotonic increase
    test=(fx_1>fx_0);

elseif test_type==1

    % Armijo sufficient increase condition
    test=(fx_1>=fx_0+spin_system.control.ls_c1*alpha*(gfx_0'*dir));

elseif test_type==2

    % Strong Wolfe curvature condition
    test=(abs(gfx_1'*dir)<=spin_system.control.ls_c2*abs(gfx_0'*dir));

elseif test_type==3

    % Positive derivative in 
    % descent direction
    test=(gfx_1'*dir>0);

end

end

% Consistency enforcement
function grumble(test_type,alpha,fx_0,fx_1,gfx_0,gfx_1,dir)
if (~isnumeric(test_type))||(~isreal(test_type))||(~isscalar(test_type))||...
   (~ismember(test_type,[0 1 2 3]))
    error('test_type must be 0, 1, 2, or 3.');
end
if ismember(test_type,[0 1])
    if isempty(fx_0)||(~isnumeric(fx_0))||(~isreal(fx_0))||(~isscalar(fx_0))
        error('fx_0 must be a real scalar for test types 0 and 1.');
    end
    if isempty(fx_1)||(~isnumeric(fx_1))||(~isreal(fx_1))||(~isscalar(fx_1))
        error('fx_1 must be a real scalar for test types 0 and 1.');
    end
end
if test_type==1
    if isempty(alpha)||(~isnumeric(alpha))||(~isreal(alpha))||(~isscalar(alpha))
        error('alpha must be a real scalar for test type 1.');
    end
end
if ismember(test_type,[1 2])
    if isempty(gfx_0)||(~isnumeric(gfx_0))||(~isreal(gfx_0))||(~iscolumn(gfx_0))
        error('gfx_0 must be a real column vector for test types 1 and 2.');
    end
end
if ismember(test_type,[2 3])
    if isempty(gfx_1)||(~isnumeric(gfx_1))||(~isreal(gfx_1))||(~iscolumn(gfx_1))
        error('gfx_1 must be a real column vector for test types 2 and 3.');
    end
end
if ismember(test_type,[1 2 3])
    if isempty(dir)||(~isnumeric(dir))||(~isreal(dir))||(~iscolumn(dir))
        error('dir must be a real column vector for test types 1, 2, and 3.');
    end
end
if (ismember(test_type,[1 2])&&(~isequal(size(gfx_0),size(dir))))||...
   (ismember(test_type,[2 3])&&(~isequal(size(gfx_1),size(dir))))
    error('dir and gradient vectors must have matching dimensions.');
end
end

% These are illusions of popular history which 
% a successful religion must promote: evil men
% never prosper; only the brave deserve the fair;
% honesty is the best policy; actions speak lou-
% der than words; virtue always triumphs; a good 
% deed is its own reward; any bad human can be 
% reformed; religious talismans protect one from
% demon possession; only females understand the
% ancient mysteries; the rich are doomed to un-
% happiness...
%
% Frank Herbert, in the Dune series

