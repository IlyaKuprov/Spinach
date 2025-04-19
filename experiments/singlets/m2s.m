% M2S sequence of Pileio and Levitt. Syntax:
%
%           rho=m2s(spin_system,L,Hx,Hy,rho,J,delta_v)
%
% Parameters:
%
%        L    - background Liouvillian
%
%       Hx    - X spin operator
%
%       Hy    - Y spin operator
%
%      rho    - initial state vector
%      
%        J    - J-coupling (Hz)
%
%  delta_v    - Zeeman frequency difference (Hz)
%
% Outputs:
%
%      rho    - final state vector.
%
% a.j.allami@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=m2s.m>

function rho=m2s(spin_system,L,Hx,Hy,rho,J,delta_v)

% Check consistency
grumble(L,Hx,Hy,rho,J,delta_v);

% Evolution time
t=1/(4*sqrt(J^2+delta_v^2));

% Repetition count
m1=floor(pi*J/(2*delta_v));
if mod(m1,2)~=0, m1=m1+1; end

% Pulse sequence
rho=step(spin_system,Hy,rho,pi/2);
for n=1:m1
    rho=step(spin_system,L,rho,t);
    rho=step(spin_system,Hx,rho,pi);
    rho=step(spin_system,L,rho,t);
end
rho=step(spin_system,Hx,rho,pi/2);
rho=step(spin_system,L,rho,t);
for n=1:m1/2
    rho=step(spin_system,L,rho,t);
    rho=step(spin_system,Hx,rho,pi);
    rho=step(spin_system,L,rho,t);
end

end

% Consistency enforcement
function grumble(L,Hx,Hy,rho,J,delta_v)
if (~isnumeric(L))||(~isnumeric(Hx))||(~isnumeric(Hy))||...
   (~ismatrix(L))||(~ismatrix(Hx))||(~ismatrix(Hy))
    error('L, Hx and Hy arguments must be matrices.');
end
if (~all(size(L)==size(Hx)))||(~all(size(Hx)==size(Hy)))
    error('L, Hx and Hy matrices must have the same dimension.');
end
if size(rho,1)~=size(L,2)
    error('dimensions of rho and L must be consistent.');
end
if (~isnumeric(J))||(~isreal(J))||(~isscalar(J))
    error('J must be a real scalar.');
end
if (~isnumeric(delta_v))||(~isreal(delta_v))||(~isscalar(delta_v))
    error('delta_v must be a real scalar.');
end
end

% Morality is herd instinct in the individual.
%
% Friedrich Nietzsche

