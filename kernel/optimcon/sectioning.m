% Docs header here.

function [alpha,fx_1,gfx_1,exitflag,data] = sectioning(cost_function,A,B,x_0,fx_0,...
                                                       gfx_0,dir,data,spin_system)

% Frost incoming gradient and direction
if ~isempty(spin_system.control.freeze)
    dir=dir.*(~spin_system.control.freeze(:));
    gfx_0=gfx_0.*(~spin_system.control.freeze(:));
end

% Enter sectioning loop
while true
    
    % Pick alpha in reduced bracket
    End_A = A.alpha + spin_system.control.ls_tau2*(B.alpha-A.alpha);
    End_B = B.alpha - spin_system.control.ls_tau3*(B.alpha-A.alpha);
    
    % Minimize cubic interpolant
    alpha = cubic_interp(End_A,End_B,A.alpha,B.alpha,A.fx,A.gfx'*dir,B.fx,B.gfx'*dir);
    
    % No acceptable point could be found
    if (abs( (alpha - A.alpha)*(A.gfx'*dir) ) <= eps(max(1,abs(fx_0))))
        exitflag = -2; return;
    end
    
    % Calculate value and gradient of current alpha
    [data,fx_1,gfx_1]=objeval(x_0+alpha*dir,cost_function,data,spin_system);

    % Frost the new gradient
    if ~isempty(spin_system.control.freeze)
        gfx_1=gfx_1.*(~spin_system.control.freeze(:));
    end

    % Store current bracket position of A
    Tmp=A;
    
    % Update the current brackets
    if (~alpha_conds(1,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system))||...
       (~alpha_conds(0,alpha,A.fx,fx_1,A.gfx,gfx_1,dir,spin_system))
        
        % Update bracket B to current alpha
        B.alpha = alpha; B.fx = fx_1; B.gfx = gfx_1;
        
    else
        
        % Wolfe conditions, if true then acceptable point found
        if alpha_conds(2,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system)
            exitflag = 0; return;
        end
        
        % Update bracket A
        A.alpha = alpha; A.fx = fx_1;  A.gfx = gfx_1;
        
        % B becomes old bracket A;
        if (A.alpha - B.alpha)*(gfx_1'*dir) >= 0, B=Tmp; end
        
    end
    
    % No acceptable point could be found
    if abs((alpha-A.alpha)*(A.gfx'*dir))<sqrt(eps)
        alpha=A.alpha; fx_1=A.fx;
        gfx_1=A.gfx; exitflag=-2; return;
    end
    
end

end

% A gentleman a few rows in front of us took grave
% exception to the behaviour of an opposing player
% and identified him, very loudly, as the author 
% of the Critique of Pure Reason - repeatedly and
% with venom.
%
% Rod Liddle, in The Spectator

