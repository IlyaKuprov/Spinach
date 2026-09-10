% Wigner function of a bosonic mode state given as a density matrix
% in a truncated Fock basis, evaluated at the specified points of the
% phase space through the displaced parity operator (Royer, Phys.
% Rev. A 15, 449, 1977):
%
%          W(alpha)=(2/pi)*trace(rho*D(alpha)*P*D(alpha)')
%
% where D(alpha)=expm(alpha*a'-conj(alpha)*a) is the displacement
% operator and P=expm(1i*pi*a'*a) is the photon number parity ope-
% rator, both built from the ladder operators truncated to the di-
% mension of the density matrix. Syntax:
%
%                      W=wigner_fock(rho,alpha)
%
% Parameters:
%
%    rho    - density matrix of the mode in the Fock basis
%             with the levels in ascending order, [n x n]
%
%    alpha  - complex phase space coordinates, an array of
%             any size; the real part is the position qua-
%             drature and the imaginary part is the momen-
%             tum quadrature in the units where a coherent
%             state |beta> has W=(2/pi)*exp(-2*|alpha-beta|^2)
%
% Outputs:
%
%    W      - Wigner function values, a real array of the
%             same size as alpha, normalised to a unit in-
%             tegral over the complex plane
%
% Note: the Fock basis must be large enough for the displaced
%       states to fit, meaning n well above |alpha|^2 at every
%       point; pad the density matrix with zero rows and co-
%       lumns when the state itself lives in a smaller space.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=wigner_fock.m>

function W=wigner_fock(rho,alpha)

% Check consistency
grumble(rho,alpha);

% Annihilation and parity operators in the truncated Fock basis
nlevels=size(rho,1); an_op=diag(sqrt(1:(nlevels-1)),1);
parity=diag((-1).^(0:(nlevels-1)));

% Displaced parity expectation values at the grid points
W=zeros(size(alpha));
for n=1:numel(alpha)
    disp_op=expm(alpha(n)*an_op'-conj(alpha(n))*an_op);
    W(n)=(2/pi)*real(trace(rho*(disp_op*parity*disp_op')));
end

end

% Consistency enforcement
function grumble(rho,alpha)
if (~isnumeric(rho))||(~ismatrix(rho))||(size(rho,1)~=size(rho,2))||(size(rho,1)<2)
    error('rho must be a square matrix of dimension at least 2.');
end
if ~isnumeric(alpha)
    error('alpha must be a numeric array.');
end
end

% Nothing is more practical than a good theory.
%
% Kurt Lewin

