% Converts a single-spin Hamiltonian back into the 
% Zeeman and quadrupolar interaction parameters that
% had been used to generate it. Syntax:
%
%                [omega,Q]=ham2nqi(H)
%
% Parameters:
%
%    H - single-spin Hamiltonian written in 
%        the Zeeman basis for a spin of any
%        multiplicity
%
% Outputs:
%
%    omega - Larmor frequencies, rad/s
%
%        Q - symmetric traceless quadrupolar 
%            coupling tensor, rad/s
%
% The outputs are returned such that:
%
%    H = omega(1)*Sx + omega(2)*Sy + omega(3)*Sz +
%             + [Sx Sy Sz]*Q*[Sx Sy Sz].';
%
% An error is produced if the Hamilonian contains
% any terms (for example, cubic) beyond those, or
% if it is not Hermitian and traceless.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ham2nqi.m>

function [omega,Q]=ham2nqi(H)

% Check consistency
grumble(H);

% Get multiplicity and basis operators
mult=size(H,1); T=irr_sph_ten(mult); S=pauli(mult);

% Get Larmor frequencies (rad/s)
omega(1)=trace(S.x'*H)/norm(S.x,'fro')^2;
omega(2)=trace(S.y'*H)/norm(S.y,'fro')^2;
omega(3)=trace(S.z'*H)/norm(S.z,'fro')^2;
omega=real(omega);

% Spin > 1/2
if mult>2

    % Get quadratic expansion coefficients
    rank2(1)=trace(T{5}'*H)/norm(T{5},'fro')^2;
    rank2(2)=trace(T{6}'*H)/norm(T{6},'fro')^2;
    rank2(3)=trace(T{7}'*H)/norm(T{7},'fro')^2;
    rank2(4)=trace(T{8}'*H)/norm(T{8},'fro')^2;
    rank2(5)=trace(T{9}'*H)/norm(T{9},'fro')^2;

    % Translate quadratic part (rad/s)
    Q=sphten2mat([],[],rank2); Q=real(Q);

else
    
    % Spin 1/2
    Q=zeros(3,3);

end

% Run an explicit reconstruction
HR=omega(1)*S.x+omega(2)*S.y+omega(3)*S.z+       ...
   Q(1,1)*S.x*S.x+Q(2,1)*S.y*S.x+Q(3,1)*S.z*S.x+ ...
   Q(1,2)*S.x*S.y+Q(2,2)*S.y*S.y+Q(3,2)*S.z*S.y+ ...
   Q(1,3)*S.x*S.z+Q(2,3)*S.y*S.z+Q(3,3)*S.z*S.z;

% Double-check reconstruction
if norm(H-HR,2)>eps()*norm(H,2)
    error('this Hamiltonian has terms of cubic or higher order.');
end

end

% Consistency enforcement
function grumble(H)
if (~isnumeric(H))||(~ishermitian(H))||...
   (abs(trace(H))>eps()*norm(H,2))||(numel(H)<4)
    error('H must be traceless, Hermitian, and at least 2x2.');
end
end

% I made a mistake in a pulse program by typing "(360) 135135" instead 
% of "(360) 135 135". I felt very guilty because I thought we'd have to
% repeat experiments. Luckily the mistake turned out to be a non-issue
% because 135135 mod 360 = 135. Malcolm asked "what god do you pray to?"
%
% Mohamed Sabba

