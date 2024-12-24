% Analytical approximation to a spin locking process. This function oblite-
% rates all spin-spin correlations and all magnetization components other
% than those along the indicated direction. Syntax:
%
%               rho=spinlock(spin_system,Lx,Ly,rho,direction)
%
% Parameters:
%
%      Lx         - X magnetization operator on the spins that
%                   should be locked
%
%      Ly         - Y magnetization operator on the spins that
%                   should be locked
%
%      rho        - state vector or a bookshelf stack thereof
%
%      direction  - direction in which the spins should be lo-
%                   cked, 'X' or 'Y'.
%
% Note: this is an approximation to what happens during a real spin locking
%       process. If you need a very accurate simulation, you would need to
%       model the spin locking explicity by adding RF terms to the system
%       Hamiltonian.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=spinlock.m>

function rho=spinlock(spin_system,Lx,Ly,rho,direction)

% Check consistency
grumble(Lx,Ly,rho,direction);

% Decide the direction
switch direction
    
    case 'X'
        
        % Destroy everything except for X magnetization
        rho=step(spin_system,Ly,rho,pi/2);
        rho=homospoil(spin_system,rho,'destroy');
        rho=step(spin_system,Ly,rho,-pi/2);
        
    case 'Y'

        % Destroy everything except for Y magnetization
        rho=step(spin_system,Lx,rho,pi/2);
        rho=homospoil(spin_system,rho,'destroy');
        rho=step(spin_system,Lx,rho,-pi/2);
        
    otherwise
        
        % Complain and bomb out
        error('unrecognized spin locking direction.');
        
end

end

% Consistency enforcement
function grumble(Lx,Ly,rho,direction)
if ~ismember(direction,{'X','Y'})
    error('direction argument can be ''X'' or ''Y''');
end
if (~isnumeric(Lx))||(size(Lx,1)~=size(Lx,2))||(~ishermitian(Lx))||...
   (~isnumeric(Ly))||(size(Ly,1)~=size(Ly,2))||(~ishermitian(Ly))
    error('Lx and Ly must be Hermitian square matrices');
end
if ~isnumeric(rho)
    error('rho parameter must be numeric');
end
end

% To every dog - his bone and cage,
% To every wolf - his teeth and rage.
%
% Victor Tsoy

