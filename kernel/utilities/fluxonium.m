% Fluxonium Hamiltonian in the truncated basis of the harmonic
% oscillator formed by the charging and the inductive energies
% of the circuit (Y. Lu, Optimal Control and Coherence Engine-
% ering for Superconducting Qubits, PhD thesis, Northwestern
% University, 2026, Eq. 2.30):
%
%    H=4*ec*n^2-ej*cos(phi-phi_e)+(el/2)*phi^2,   [phi,n]=1i
%
% where ec is the charging energy, ej is the Josephson energy,
% el=(Phi_0/2*pi)^2/L is the inductive energy of the shunt in-
% ductance L (Eq. 2.31, the factor 1/2 is explicit in the po-
% tential and is not absorbed into el), and phi_e=2*pi*Phi_e/
% Phi_0 is the reduced external flux (Eq. 2.32). Eq. 2.30 has
% the inductive term as (el/2)*(phi+phi_e)^2 and the Josephson
% term as -ej*cos(phi); the form above is its image under the
% translation of the phase origin to the minimum of the induc-
% tive potential, which does not change the spectrum. The off-
% set charge of Eq. 2.30 is dropped because the unbounded phase
% variable makes it removable by a gauge transformation. Phase
% and charge are built from the ladder operators of the linear
% oscillator 4*ec*n^2+(el/2)*phi^2 with the plasma frequency
% sqrt(8*ec*el) as
%
%          phi=phi_zpf*(b'+b),   n=1i*n_zpf*(b'-b)
%
%     phi_zpf=(2*ec/el)^(1/4),   n_zpf=(el/(32*ec))^(1/4)
%
% and the cosine is computed as a matrix function of the Her-
% mitian phase operator using the matrix exponential.
%
% Syntax:
%
%      [H,n_op,phi_op]=fluxonium(ec,ej,el,phi_e,nlevels)
%
% Parameters:
%
%   ec      - charging energy in Hz (energy over the
%             Planck constant), a positive real number
%
%   ej      - Josephson energy in Hz (energy over the
%             Planck constant), a positive real number
%
%   el      - inductive energy in Hz (energy over the
%             Planck constant), a positive real number
%
%   phi_e   - reduced external flux 2*pi*Phi_e/Phi_0
%             in radians, a real number
%
%   nlevels - number of oscillator basis states, a
%             positive integer
%
% Outputs:
%
%   H       - fluxonium Hamiltonian in rad/s (2*pi times
%             the energy in Hz), a real symmetric matrix
%             of dimension nlevels
%
%   n_op    - charge operator (Cooper pair number) in the
%             oscillator basis, a Hermitian matrix of di-
%             mension nlevels
%
%   phi_op  - phase operator in radians in the oscillator
%             basis, a real symmetric matrix of dimension
%             nlevels
%
% Note: the oscillator basis must be large enough for the lowest
%       eigenstates to be converged; nlevels of the order of 30
%       to 60 is sufficient for ej/el=5 and ec/el=1, larger ej/el
%       ratios need more states. Check the convergence by repeat-
%       ing the calculation with a bigger nlevels.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=fluxonium.m>

function [H,n_op,phi_op]=fluxonium(ec,ej,el,phi_e,nlevels)

% Check consistency
grumble(ec,ej,el,phi_e,nlevels);

% Ladder operators of the linear oscillator
A=weyl(nlevels);

% Zero-point amplitudes of the phase and the charge
phi_zpf=(2*ec/el)^(1/4); n_zpf=(el/(32*ec))^(1/4);

% Phase and charge operators
phi_op=full(phi_zpf*(A.c+A.a)); n_op=full(1i*n_zpf*(A.c-A.a));

% Cosine of the flux-shifted phase, real part drops round-off
U=expm(1i*(phi_op-phi_e*eye(nlevels))); cos_op=real(U+U')/2;

% Hamiltonian in rad/s
H=2*pi*(4*ec*n_op^2-ej*cos_op+(el/2)*phi_op^2);

end

% Consistency enforcement
function grumble(ec,ej,el,phi_e,nlevels)
if (~isnumeric(ec))||(~isreal(ec))||(~isscalar(ec))||(~isfinite(ec))||(ec<=0)
    error('ec must be a positive real number.');
end
if (~isnumeric(ej))||(~isreal(ej))||(~isscalar(ej))||(~isfinite(ej))||(ej<=0)
    error('ej must be a positive real number.');
end
if (~isnumeric(el))||(~isreal(el))||(~isscalar(el))||(~isfinite(el))||(el<=0)
    error('el must be a positive real number.');
end
if (~isnumeric(phi_e))||(~isreal(phi_e))||(~isscalar(phi_e))||(~isfinite(phi_e))
    error('phi_e must be a real number.');
end
if (~isnumeric(nlevels))||(~isreal(nlevels))||(~isscalar(nlevels))||...
   (mod(nlevels,1)~=0)||(nlevels<1)
    error('nlevels must be a positive real integer.');
end
end

% The first draft of anything is shit.
%
% Ernest Hemingway

