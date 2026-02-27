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
%     method - 'PWCL', 'PWCM', 'RKMK4', 'RKMK-DP5',
%              'RKMK-DP8', or 'LG4'; the latter one
%               has a good balance of efficiency and
%               numerical accuracy
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

function rho_b=iserstep(spin_system,LTM,rho_a,dt)

% Check consistency
grumble(LTM,rho_a,dt);

% Dereference input pointers
L=LTM{1}; t=LTM{2}; method=LTM{3}; 

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
    
    case 'LG4A' % https://doi.org/10.1088/0305-4470/39/19/S07
        
        % Match the paper notation
        A_func = @(t,rho)(-1i*L(t,rho));
        h=dt; t_n=t; Y_n=rho_a;

        % Stage 1
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
        rho_mid=step(spin_system,LMA+1i*dt*(1/6)*comm(LL,LMA),rho_a,dt/2);
        LMB=L(t+0.5*dt,rho_mid);
        
        % Estimate right point state and generator
        rho_right=step(spin_system,LMB+1i*dt*(1/6)*comm(LL,LMB),rho_a,dt);
        LR=L(t+dt,rho_right);
        
        % Get the average generator and the commutator correction
        LA=(1/6)*(LL+2*LMA+2*LMB+LR); LC=comm(LL-LR,LMA+LMB)+comm(LL,LR);

        % Take the step under RKMK4 generator
        rho_b=step(spin_system,LA+1i*dt*(1/36)*LC,rho_a,dt);

    case 'RKMK-DP5' % Fifth order Dormand-Prince RKMK

        % Take the step under DP5 tableau
        rho_b=RKMK('DP5',spin_system,L,t,rho_a,dt);

    case 'RKMK-DP8' % Eighth order Dormand-Prince RKMK

        % Take the step under DP8 tableau
        rho_b=RKMK('DP8',spin_system,L,t,rho_a,dt);

    case 'RKMK-RKF45' % Fifth order Fehlberg RKMK

        % Take the step under RKF45 tableau
        rho_b=RKMK('RKF45',spin_system,L,t,rho_a,dt);

    otherwise
        
        % Complain and bomb out
        error('unknown propagation method');
        
end

end

% Consistency enforcement
function grumble(LTM,rho_a,dt)
if ~iscell(LTM)
    error('second input must be a cell array.');
end
if ~isa(LTM{1},'function_handle')
    error('evolution generator must be a function handle.');
end
if (~isnumeric(LTM{2}))||(~isscalar(LTM{2}))
    error('current time must be a numeric scalar.');
end
if ~ischar(LTM{3})
   error('method must be a character string.')
end
if (~isnumeric(dt))||(~isnumeric(rho_a))
    error('rho_a and dt must be numeric.');
end
if ~isscalar(dt)
    error('dt must be a scalar.');
end
end

% Explicit Runge-Kutta-Munthe-Kaas step
function rho_b=RKMK(order_or_name,spin_system,L,t,rho_a,dt,varargin)

% Process built-in method names
if ischar(order_or_name)||isstring(order_or_name)
    method_name=char(order_or_name);
    [A,b,c,q]=rk_tableau(method_name);
    q_dexp=q-1;
else
    q=order_or_name;
    A=[]; b=[]; c=[];
    q_dexp=q-1;
end

% Process optional arguments
k=1;
while k<=numel(varargin)
    key=varargin{k};
    val=varargin{k+1};
    switch lower(key)
        case 'a'
            A=val;
        case 'b'
            b=val;
        case 'c'
            c=val;
        case 'q_dexp'
            q_dexp=val;
        otherwise
            error('RKMK:BadOption','unknown option "%s".',key);
    end
    k=k+2;
end

% Check tableau presence
if isempty(A)||isempty(b)||isempty(c)
    if ~(ischar(order_or_name)||isstring(order_or_name))
        error('RKMK:TableauMissing','numeric order requires A,b,c.');
    end
end

% Get the number of stages
n_stages=numel(b);

% Build scaling factor
prop_scale=1i/dt;

% Get commutator function
comm_fun=@(X,Y)X*Y-Y*X;

% Get Bernoulli coefficients
bern_coeffs=bernoulli_B(q_dexp);

% Preallocate stage generators
K=cell(n_stages,1);

% Run stage loop
for n=1:n_stages

    % Build Lie algebra element for this stage
    u_n=0;
    if n>1
        for m=1:n-1
            if A(n,m)~=0
                u_n=u_n+dt*A(n,m)*K{m};
            end
        end
    end

    % Propagate to stage state
    if isequal(u_n,0)
        rho_n=rho_a;
    else
        rho_n=step(spin_system,prop_scale*u_n,rho_a,dt);
    end

    % Evaluate algebra field at stage node
    k_raw=-1i*L(t+c(n)*dt,rho_n);

    % Apply truncated dexp^{-1} map
    if isequal(u_n,0)||q_dexp==0
        K{n}=k_raw;
    else
        K{n}=dexpinv(u_n,k_raw,q_dexp,comm_fun,bern_coeffs);
    end

end

% Accumulate final algebra element
v=0;
for n=1:n_stages
    if b(n)~=0
        v=v+dt*b(n)*K{n};
    end
end

% Propagate final state
if isequal(v,0)
    rho_b=rho_a;
else
    rho_b=step(spin_system,prop_scale*v,rho_a,dt);
end

end

% Truncated Bernoulli series for dexp^{-1}
function w=dexpinv(u,v,q_dexp,comm_fun,bern_coeffs)

% Initialise the series and first adjoint
w=v;
ad_term=v;

% Build the truncated series
for n=1:q_dexp
    ad_term=comm_fun(u,ad_term);
    Bn=bern_coeffs(n+1);
    if Bn~=0
        w=w+(Bn/factorial(n))*ad_term;
    end
end

end

% Bernoulli numbers B_0..B_n for dexp^{-1}
function B=bernoulli_B(n)

% Set available Bernoulli numbers
B_all=[1,-1/2,1/6,0,-1/30,0,1/42,0,-1/30,0,5/66];

% Make sure the requested order is available
if n>10
    error('RKMK:Bernoulli','extend bernoulli_B for n=%d.',n);
end

% Return requested prefix
B=B_all(1:n+1);

end

% Explicit RK tableaux for named methods
function [A,b,c,q]=rk_tableau(name)

% Standardise method name
name=upper(strtrim(name));

% Select tableau
switch name

    case {'DP5','DORMANDPRINCE54','DOPRI5'}

        % Set global order
        q=5;

        % Set stage nodes
        c=[0;1/5;3/10;4/5;8/9;1;1];

        % Set stage matrix
        A=zeros(7,7);
        A(2,1)=1/5;
        A(3,1)=3/40;       A(3,2)=9/40;
        A(4,1)=44/45;      A(4,2)=-56/15;        A(4,3)=32/9;
        A(5,1)=19372/6561; A(5,2)=-25360/2187;   A(5,3)=64448/6561; A(5,4)=-212/729;
        A(6,1)=9017/3168;  A(6,2)=-355/33;       A(6,3)=46732/5247; A(6,4)=49/176;   A(6,5)=-5103/18656;
        A(7,1)=35/384;     A(7,3)=500/1113;      A(7,4)=125/192;    A(7,5)=-2187/6784; A(7,6)=11/84;

        % Set output weights
        b=[35/384;0;500/1113;125/192;-2187/6784;11/84;0];

    case {'RKF45','FEHLBERG45','FEHLBERG4(5)'}

        % Set global order
        q=5;

        % Set stage nodes
        c=[0;1/4;3/8;12/13;1;1/2];

        % Set stage matrix
        A=zeros(6,6);
        A(2,1)=1/4;
        A(3,1)=3/32;       A(3,2)=9/32;
        A(4,1)=1932/2197;  A(4,2)=-7200/2197;    A(4,3)=7296/2197;
        A(5,1)=439/216;    A(5,2)=-8;            A(5,3)=3680/513;   A(5,4)=-845/4104;
        A(6,1)=-8/27;      A(6,2)=2;             A(6,3)=-3544/2565; A(6,4)=1859/4104; A(6,5)=-11/40;

        % Set output weights
        b=[16/135;0;6656/12825;28561/56430;-9/50;2/55];

    case {'DP8','DORMANDPRINCE87','DOPRI8'}

        % Set global order
        q=8;

        % Set stage nodes
        c=zeros(13,1);
        c(1)=0;
        c(2)=1/18;
        c(3)=1/12;
        c(4)=1/8;
        c(5)=5/16;
        c(6)=3/8;
        c(7)=59/400;
        c(8)=93/200;
        c(9)=5490023248/9719169821;
        c(10)=13/20;
        c(11)=1201146811/1299019798;
        c(12)=1;
        c(13)=1;

        % Set stage matrix
        A=zeros(13,13);
        A(2,1)=1/18;
        A(3,1)=1/48; A(3,2)=1/16;
        A(4,1)=1/32; A(4,3)=3/32;
        A(5,1)=5/16; A(5,3)=-75/64; A(5,4)=75/64;
        A(6,1)=3/80; A(6,4)=3/16; A(6,5)=3/20;
        A(7,1)=29443841/614563906;
        A(7,4)=77736538/692538347;
        A(7,5)=-28693883/1125000000;
        A(7,6)=23124283/1800000000;
        A(8,1)=16016141/946692911;
        A(8,4)=61564180/158732637;
        A(8,5)=22789713/633445777;
        A(8,6)=545815736/2771057229;
        A(8,7)=-180193667/1043307555;
        A(9,1)=39632708/573591083;
        A(9,4)=-433636366/683701615;
        A(9,5)=-421739975/2616292301;
        A(9,6)=100302831/723423059;
        A(9,7)=790204164/839813087;
        A(9,8)=800635310/3783071287;
        A(10,1)=246121993/1340847787;
        A(10,4)=-37695042795/15268766246;
        A(10,5)=-309121744/1061227803;
        A(10,6)=-12992083/490766935;
        A(10,7)=6005943493/2108947869;
        A(10,8)=393006217/1396673457;
        A(10,9)=123872331/1001029789;
        A(11,1)=-1028468189/846180014;
        A(11,4)=8478235783/508512852;
        A(11,5)=1311729495/1432422823;
        A(11,6)=-10304129995/1701304382;
        A(11,7)=-48777925059/3047939560;
        A(11,8)=15336726248/1032824649;
        A(11,9)=-45442868181/3398467696;
        A(11,10)=3065993473/597172653;
        A(12,1)=185892177/718116043;
        A(12,4)=-3185094517/667107341;
        A(12,5)=-477755414/1098053517;
        A(12,6)=-703635378/230739211;
        A(12,7)=5731566787/1027545527;
        A(12,8)=5232866602/850066563;
        A(12,9)=-4093664535/808688257;
        A(12,10)=3962137247/1805957418;
        A(12,11)=65686358/487910083;
        A(13,1)=403863854/491063109;
        A(13,4)=-5068492393/434740067;
        A(13,5)=-411421997/543043805;
        A(13,6)=652783627/914296604;
        A(13,7)=11173962825/925320556;
        A(13,8)=-13158990841/6184727034;
        A(13,9)=3936647629/1978049680;
        A(13,10)=-160528059/685178525;
        A(13,11)=248638103/1413531060;
        A(13,12)=0;

        % Set output weights
        b=zeros(13,1);
        b(1)=14005451/335480064;
        b(2)=0;
        b(3)=0;
        b(4)=0;
        b(5)=0;
        b(6)=-59238493/1068277825;
        b(7)=181606767/758867731;
        b(8)=561292985/797845732;
        b(9)=-1041891430/1371343529;
        b(10)=760417239/1151165299;
        b(11)=118820643/751138087;
        b(12)=-528747749/2220607170;
        b(13)=1/4;

    otherwise

        % Complain and bomb out
        error('rk_tableau:UnknownMethod','unknown method "%s".',name);

end

end

% The last bunch of pickets were carrying signs that 
% said "Make Love, Not War". The only trouble was 
% they didn't look capable of doing either.
%
% Ronald Reagan

