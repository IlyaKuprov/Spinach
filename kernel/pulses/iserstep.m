% Lie-group and Runge-Kutta-Munthe-Kaas solvers for the Lie equa-
% tion. LG methods are implementations of Equation A.1, with mi-
% nor typos fixed, from 
%
%          http://dx.doi.org/10.1088/0305-4470/39/19/S07
%
% The key difference from step() function is that the Liouvillian
% can depend on the density matrix. Syntax:
%
%         rho_b=iserstep(spin_system,{L,t,method},rho_a,dt)
%
% Parameters:
%
%     spin_system - Spinach data structure from create.m
%                   and basis.m constructors
%
%     L - a handle to a function L(t,rho) that must take 
%         time and state vector, and return the evolution
%         generator (in rad/s) of the Lie equation:
%
%                   d_rho/d_t = -i*L(t,rho)*rho
%
%     rho_a - state vector at the start of the evolution 
%             period
%
%     t   - time at the start of the evolution, seconds
%
%     dt  - evolution time step, seconds
%
%     method - 'PWCL', 'PWCM', 'RKMK4', or 'LG4',
%              the latter one is recommended
%
% Outputs:
%
%     rho_b - state vector at the end of the evolution 
%             time step
%
% ilya.kuprov@weizmann.ac.il
% a.graham@soton.ac.uk
% a.acharya@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=iserstep.m>

function rho_b=iserstep(spin_system,L_t_method,rho_a,dt)
L=L_t_method{1};
t=L_t_method{2};
method=L_t_method{3}; 

% Check consistency
grumble(L,rho_a,t,dt,method);

% Decide method
switch method
    
    case 'PWCL'  % Piecewise-constant method

        % Left generator from left state
        LL=L(t,rho_a);
        
        % Assume the generator stays constant
        rho_b=step(spin_system,LL,rho_a,dt);

    case 'LG2'   % Second order Lie-group method

        % Left generator from left state
        LL=L(t,rho_a);
        
        % Estimate midpoint state and generator
        rho_mid=step(spin_system,LL,rho_a,dt/2);
        LM=L(t+0.5*dt,rho_mid);

        % Assume the generator stays constant
        rho_b=step(spin_system,LM,rho_a,dt);

    case 'LG4'   % Fourth order Lie-group method

        % Left generator from left state
        LL=L(t,rho_a);

        % Estimate midpoint state and generator
        rho_mid=step(spin_system,LL,rho_a,dt/2);
        LM=L(t+0.5*dt,rho_mid);
        
        % Estimate endpoint using midpoint generator
        rho_b=step(spin_system,LM,rho_a,dt);
        
        % Take the step using fourth-order Lie method
        rho_b=step(spin_system,{LL,LM,L(t+dt,rho_b)},rho_a,dt);
    
    case 'LG4A'
        A_func = @(t,rho) -1i*L(t,rho);
        h=dt; t_n=t; Y_n=rho_a;
        % Stage 1
        u1 = 0;
        k1 = h * A_func(t_n, Y_n);
        Q1 = k1;

        % Stage 2
        u2 = 0.5 * Q1;
        k2 = h * A_func(t_n + h/2, expm(u2) * Y_n);
        Q2 = k2 - k1;

        % Stage 3
        u3 = 0.5 * Q1 + 0.25 * Q2;
        k3 = h * A_func(t_n + h/2, expm(u3) * Y_n);
        Q3 = k3 - k2;

        % Stage 4
        u4 = Q1 + Q2;
        k4 = h * A_func(t_n + h, expm(u4) * Y_n);
        Q4 = k4 - 2*k2 + k1;

        % Stage 5 (Note: [A, B] is the commutator A*B - B*A)
        u5 = 0.5*Q1 + 0.25*Q2 + (1/3)*Q3 - (1/24)*Q4 - (1/48)*(Q1*Q2 - Q2*Q1);
        k5 = h * A_func(t_n + h/2, expm(u5) * Y_n);
        Q5 = k5 - k2;

        % Stage 6
        u6 = Q1 + Q2 + (2/3)*Q3 + (1/6)*Q4 - (1/6)*(Q1*Q2 - Q2*Q1);
        k6 = h * A_func(t_n + h, expm(u6) * Y_n);
        Q6 = k6 - 2*k2 + k1;

        % Final exponent v
        v = Q1 + Q2 + (2/3)*Q5 + (1/6)*Q6 - (1/6)*(Q1*Q2 - (Q2 - Q3 + Q5 + 0.5*Q6)*Q1);

        % Step update
        Y_next = expm(v) * Y_n;
        rho_b = Y_next;

    case 'RKMK4' % Fourth order RKMK method

        % Left generator from left state
        LL=L(t,rho_a);
            
        % Estimate midpoint state and generator
        rho_mid=step(spin_system,LL,rho_a,dt/2);
        LMA=L(t+0.5*dt,rho_mid);
        
        % Re-estimate midpoint state and generator
        rho_mid=step(spin_system,LMA+1i*dt*(1/6)*(LL*LMA-LMA*LL),rho_a,dt/2);
        LMB=L(t+0.5*dt,rho_mid);
        
        % Estimate right point state and generator
        rho_right=step(spin_system,LMB+1i*dt*(1/6)*(LL*LMB-LMB*LL),rho_a,dt);
        LR=L(t+dt,rho_right);
        
        % Get the average generator and the commutator correction
        LA=(1/6)*(LL+2*LMA+2*LMB+LR); LC=comm(LL-LR,LMA+LMB)+comm(LL,LR);

        % Take the step under RKMK4 generator
        rho_b=step(spin_system,LA+1i*dt*(1/36)*LC,rho_a,dt);

    otherwise
        
        % Complain and bomb out
        error('unknown propagation method');
        
end

end

% Consistency enforcement
function grumble(L,rho_a,t,dt,method)
if ~isa(L,'function_handle')
    error('L must be a function handle.');
end
if (~isnumeric(t))||(~isnumeric(dt))||(~isnumeric(rho_a))
    error('rho_a, t, and dt must be numeric.');
end
if (~isscalar(t))||(~isscalar(dt))
    error('t and dt must be scalars.');
end
if ~ischar(method)
   error('method must be a character string.')
end
end

% The last bunch of pickets were carrying signs that 
% said "Make Love, Not War". The only trouble was 
% they didn't look capable of doing either.
%
% Ronald Reagan

