% Docs header here.

function [A,B,alpha,fx,gfx,next_act,data]=bracketing(cost_function,alpha,dir,x_0,...
                                                     fx_0,gfx_0,data,spin_system)

% Take a look at the proposed new point
[data,fx_2,gfx_2]=objeval(x_0+dir,cost_function,data,spin_system);

% Initialise brackets A and B
A.alpha=[]; A.fx=[]; A.gfx=[];
B.alpha=[]; B.fx=[]; B.gfx=[];

% Frost incoming gradients and direction
if ~isempty(spin_system.control.freeze)
    dir=dir.*(~spin_system.control.freeze(:));
    gfx_0=gfx_0.*(~spin_system.control.freeze(:));
    gfx_2=gfx_2.*(~spin_system.control.freeze(:));
end

% Get the history going
fx=fx_0;   fx_1=fx_0;
gfx=gfx_0; gfx_1=gfx_0;
alpha_1=0; alpha_2=alpha;

% Enter the bracketing loop
while true
    
    % Bracket located - case 1 (Wolfe conditions)
    if (~alpha_conds(1,alpha_2,fx_0,fx_2,gfx_0,[],dir,spin_system))||...
       (~alpha_conds(0,[],fx_1,fx_2,[],[],[],spin_system))
   
        % Set the bracket values
        A.alpha=alpha_1; A.fx=fx_1; A.gfx=gfx_1;
        B.alpha=alpha_2; B.fx=fx_2; B.gfx=gfx_2;
        
        % Proceed to sectioning
        next_act='sectioning'; return;
        
    end
    
    % Acceptable step length found
    if alpha_conds(2,[],[],[],gfx_0,gfx_2,dir,spin_system)
        
        % Store the found alpha values
        alpha=alpha_2; fx=fx_2; gfx=gfx_2;
        
        % No sectioning required
        next_act='none';  return;
        
    end
    
    % Bracket located - case 2
    if ~alpha_conds(3,[],[],[],[],gfx_2,dir,spin_system)
        
        % Set the bracket values
        A.alpha=alpha_2; A.fx=fx_2; A.gfx=gfx_2;
        B.alpha=alpha_1; B.fx=fx_1; B.gfx=gfx_1;
        
        % Finished bracketing phase
        next_act='sectioning'; return;
        
    end
    
    % Get bracket endpoints
    br_end_pt_A=2*alpha_2-alpha_1;
    br_end_pt_B=alpha_2+spin_system.control.ls_tau1*(alpha_2-alpha_1);
    
    % Maximise cubic interpolant
    alpha_new=cubic_interp(br_end_pt_A,br_end_pt_B,alpha_1,alpha_2,...
                           fx_1,gfx_1'*dir,fx_2,gfx_2'*dir);
    
    % Update history
    alpha_1=alpha_2; alpha_2=alpha_new;
    fx_1=fx_2; gfx_1=gfx_2; x_1=x_0+alpha_2*dir;
    
    % Calculate value and gradient for current alpha
    [data,fx_2,gfx_2]=objeval(x_1,cost_function,data,spin_system);

    % Frost the new gradient
    if ~isempty(spin_system.control.freeze)
        gfx_2=gfx_2.*(~spin_system.control.freeze(:));
    end
    
end

end

% Одной рукой бунтую, другой пишу донос.
%
% Михаил Щербаков

