% Performs a line search to find an appropriate optimisation step length
% in the specified direction. Based on fminlbfgs.m code from D. Kroon,
% University of Twente (Nov 2010). Syntax:
%
%  [alpha,fx_1,gfx_1,exitflag,data]=...
%          linesearch(spin_system,cost_function,d_0,x_0,fx_0,gfx_0,data)
%
% Parameters:
%
%    cost_function  - objective function handle
%
%    d_0    - proposed step vector
%
%    x_0    - current point
%
%    fx_0   - function value at the current point
%
%    gfx_0  - gradient at the current point
%
%    data   - diagnostic data structure
%
% Outputs:
%
%    alpha     - optimum step length multiplier
%
%    fx_1      - function value at the resulting point
%
%    gfx_1     - gradient at the resulting point
%
%    exitflag  - termination message
%
%    data      - updated disgnostic data structure
%
% david.goodwin@inano.au.dk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=linesearch.m>

function [alpha,fx_1,gfx_1,exitflag,data]=...
          linesearch(spin_system,cost_function,d_0,x_0,fx_0,gfx_0,data)

% Check consistency
grumble(cost_function);
      
% Step length guess
if strcmp(spin_system.control.method,'lbfgs')&&(data.count.iter==1)
    
    % Conservative first step for LBFGS
    alpha=min(1/norm(gfx_0,inf),5);
    
else
    
    % Hessian knows best
    alpha=1;
    
end

% Take a look at the proposed new point
[data,fx_1,gfx_1]=objeval(x_0+alpha*d_0,cost_function,data,spin_system);

% Find a bracket [A B] of acceptable points
[A,B,alpha,fx_1,gfx_1,exitflag,data]=bracketing(cost_function,alpha,d_0,x_0,...
                                     fx_0,gfx_0,fx_1,gfx_1,data,spin_system);

% Run sectioning
if exitflag==2
    
    % Find acceptable point within bracket    
    [alpha,fx_1,gfx_1,exitflag,data]=sectioning(cost_function,A,B,x_0,fx_0,...
                                     gfx_0,d_0,data,spin_system);
                                 
end

end

%% Bracketing subroutine
function [A,B,alpha,fx,gfx,exitflag,data] = bracketing(cost_function,alpha,d_0,x_0,fx_0,...
                                                       gfx_0,fx_2,gfx_2,data,spin_system)

% Initialise bracket A and bracket B
A.alpha = []; A.fx = []; A.gfx = [];
B.alpha = []; B.fx = []; B.gfx = [];

% Get the history going
fx=fx_0;   fx_1=fx_0;
gfx=gfx_0; gfx_1=gfx_0;
alpha_1=0; alpha_2=alpha;

% Enter the bracketing loop
while true
    
    % Bracket located - case 1 (Wolfe conditions)
    if (~alpha_conditions(1,alpha,fx_0,fx_2,gfx_0,[],d_0,spin_system))||...
       (~alpha_conditions(0,[],fx_0,fx_2,[],[],[],spin_system))
   
        % Set the bracket values
        A.alpha = alpha_1; A.fx = fx_1;  A.gfx = gfx_1;
        B.alpha = alpha_2; B.fx = fx_2;  B.gfx = gfx_2;
        
        % Finished bracketing phase
        exitflag = 2; return;
        
    end
    
    % Acceptable steplength found
    if alpha_conditions(2,alpha,fx_0,fx_2,gfx_0,gfx_2,d_0,spin_system)
        
        % Store the found alpha values
        alpha=alpha_2; fx=fx_2; gfx=gfx_2;
        
        % Finished bracketing phase, no need to call sectioning phase
        exitflag = [];  return;
        
    end
    
    % Bracket located - case 2
    if ~alpha_conditions(3,[],[],[],[],gfx_2,d_0,spin_system)
        
        % Set the bracket values
        A.alpha = alpha_2; A.fx = fx_2;  A.gfx = gfx_2;
        B.alpha = alpha_1; B.fx = fx_1;  B.gfx = gfx_1;
        
        % Finished bracketing phase
        exitflag  = 2; return;
        
    end
    
    % Get the endpoints
    brcktEndpntA = 2*alpha_2-alpha_1;
    brcktEndpntB = alpha_2+spin_system.control.ls_tau1*(alpha_2-alpha_1);
    
    % Minimize cubic interpolant
    alpha_new = cubic_interpolation(brcktEndpntA,brcktEndpntB,alpha_1,alpha_2,...
                                    fx_1,gfx_1'*d_0,fx_2,gfx_2'*d_0);
    
    % Update history
    alpha_1=alpha_2; 
    alpha_2=alpha_new;
    fx_1=fx_2; gfx_1=gfx_2;
    x_1=x_0+alpha_2*d_0;
    
    % Calculate value and gradient for current alpha
    [data,fx_2,gfx_2]=objeval(x_1,cost_function,data,spin_system);
    
end

end

%% Sectioning subroutine
function [alpha,fx_1,gfx_1,exitflag,data] = sectioning(cost_function,A,B,x_0,fx_0,...
                                                       gfx_0,d_0,data,spin_system)

% Worst-case returns are empty
fx_1=[]; gfx_1=[];

% Enter sectioning loop
while true
    
    % Pick alpha in reduced bracket
    End_A = A.alpha + min(spin_system.control.ls_tau2,...
                          spin_system.control.ls_c2)*(B.alpha - A.alpha);
    End_B = B.alpha - spin_system.control.ls_tau3*(B.alpha - A.alpha);
    
    % Minimize cubic interpolant
    alpha = cubic_interpolation(End_A,End_B,A.alpha,B.alpha,A.fx,A.gfx'*d_0,B.fx,B.gfx'*d_0);
    
    % No acceptable point could be found
    if (abs( (alpha - A.alpha)*(A.gfx'*d_0) ) <= eps(max(1,abs(fx_0))))
        exitflag = -2; return;
    end
    
    % Calculate value and gradient of current alpha
    [data,fx_1,gfx_1]=objeval(x_0+alpha*d_0,cost_function,data,spin_system);
    
    % Store current bracket position of A
    Tmp=A;
    
    % Update the current brackets
    if (~alpha_conditions(1,alpha,fx_0,fx_1,gfx_0,gfx_1,d_0,spin_system))||...
       (~alpha_conditions(0,alpha,A.fx,fx_1,A.gfx,gfx_1,d_0,spin_system))
        
        % Update bracket B to current alpha
        B.alpha = alpha; B.fx = fx_1; B.gfx = gfx_1;
        
    else
        
        % Wolfe conditions, if true then acceptable point found
        if alpha_conditions(2,alpha,fx_0,fx_1,gfx_0,gfx_1,d_0,spin_system)
            exitflag = []; return;
        end
        
        % Update bracket A
        A.alpha = alpha; A.fx = fx_1;  A.gfx = gfx_1;
        
        % B becomes old bracket A;
        if (A.alpha - B.alpha)*(gfx_1'*d_0) >= 0, B=Tmp; end
        
    end
    
    % No acceptable point could be found
    if (abs(B.alpha-A.alpha) < eps), exitflag = -2; return, end
    
end

end

function [alpha,fx]=cubic_interpolation(End_A,End_B,alpha_A,alpha_B,f_A,dir_deriv_A,f_B,dir_deriv_B)

% Get the coefficients of the cubic polynomial
c1=-2*(f_B-f_A)+(dir_deriv_A+dir_deriv_B)*(alpha_B-alpha_A);
c2= 3*(f_B-f_A)-(2*dir_deriv_A+dir_deriv_B)*(alpha_B-alpha_A);
c3=(alpha_B-alpha_A)*dir_deriv_A; c4=f_A;

% Convert bounds to the z-space
bounds = ([End_A End_B ]-alpha_A)./(alpha_B - alpha_A);

% Find minima and maxima from the roots of the derivative
sPoints = roots([3*c1 2*c2 1*c3]);

% Remove imaginary and points outside range and make vector with solutions
sPoints(imag(sPoints)~=0)=[];
sPoints(sPoints<min(bounds))=[];
sPoints(sPoints>max(bounds))=[];
sPoints=[min(bounds) sPoints(:)' max(bounds)];

% Select the global minimum point
[fx,k]=max(polyval([c1 c2 c3 c4],sPoints));

% Add the offset and scale back from [0..1] to the alpha domain
alpha = alpha_A + sPoints(k)*(alpha_B - alpha_A);

end

function test=alpha_conditions(test_type,alpha,fx_0,fx_1,gfx_0,gfx_1,d_0,spin_system)

% Select test type
if test_type==0
    
    % Gradient free test - ensure increase
    test=(fx_1 > fx_0);
    
elseif test_type==1
    
    % Armijo condition - ensure sufficient increase
    test=(fx_1 >= fx_0 + spin_system.control.ls_c1*alpha*(gfx_0'*d_0));
    
elseif test_type==2
    
    % Strong Wolfe (Curvature) condition - avoid underestimated step length
    test=(abs(gfx_1'*d_0) <= spin_system.control.ls_c2*abs(gfx_0'*d_0));
    
elseif test_type==3
    
    % Gradient test - ensure an ascent direction
    test=(gfx_1'*d_0 > 0);
    
end

end

% Consistency enforcement
function grumble(cost_function)
if ~isa(cost_function,'function_handle')
    error('cost_function must be a function handle.');
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

