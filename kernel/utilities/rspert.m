% Rayleigh-Schrodinger perturbation theory to arbitrary order, Eqs 
% 2.21-2.23 from Stefan Stoll's PhD thesis, with the typo fixed in
% the numerator of Eq 2.21. Syntax:
%
%                    [Ep,Vp]=rspert(E0,H1,order)
%
% Parameters:
%
%     E0    - eigenvalues of H0, a column vector of real 
%             numbers
%
%     H1    - perturbation, written in the basis that di-
%             agonalises H0
%
%     order - order of perturbation theory to be used, 6
%             is the sensible maximum
%
% Outputs:
%
%     Ep    - eigenvalues of H0+H1 to the specified order,
%             a vector of reals, not necessarily sorted in
%             the same way as the input
%
%     Vp    - normalised eigenvectors of H0+H1 to the spe-
%             cified order in perturbation theory, a squa-
%             re unitary matrix with eigenvectors in cols
%             in the same order as the eigenvalues in Ep
%
% Notes: there must be no degeneracies in H0; H1 must be Hermitian,
%        the theory only converges when norm(H1,2) is much smaller
%        than the smallest energy gap in H0; numerical artefacts
%        appear beyond sixth order; complexity is linear in the or-
%        der and cubic in the matrix dimension.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rspert.m>

function [Ep,Vp]=rspert(E0,H1,order)

% Check consistency
grumble(E0,H1,order);

% Auxiliary arrays
E=cell(order,1); V=cell(order,1);

% Reciprocal energy differences
Q=1./(E0'-E0); Q(logical(speye(size(Q))))=0;
if any(~isfinite(Q),'all')
    error('H0 has degenerate energy levels.');
end

% First order
E{1}=diag(H1); V{1}=Q.*H1;
    
% Higher orders
for k=2:order
    
    % Expensive
    R=H1*V{k-1};
    
    % Eigenvalues
    E{k}=real(diag(R));
    
    % Eigenvectors
    V{k}=R;
    for m=1:(k-1)
        V{k}=V{k}-V{k-m}.*E{m}';
    end
    V{k}=Q.*V{k};
    
end

% Summation
Ep=E0; for n=1:order, Ep=Ep+E{n}; end
Vp=eye(size(H1)); for n=1:order, Vp=Vp+V{n}; end 

% Normalisation
Vp=Vp./sqrt(sum(abs(Vp).^2,1));

end

% Consistency enforcement
function grumble(E0,H1,order)
if (~isnumeric(E0))||(~isreal(E0))||(~iscolumn(E0))
    error('E0 must be a real column vector.');
end
if (~isnumeric(H1))||(~ishermitian(H1))
    error('H1 must be a Hermitian matrix.');
end
if numel(E0)~=size(H1,1)
    error('dimensions of E0 and H1 must be consistent.');
end
if (~isnumeric(order))||(~isreal(order))||...
   (~isscalar(order))||(mod(order,1)~=0)||(order<1)
    error('order must be a positive real integer.');
end
end

% "How can you sleep at night when the Universe is expanding?!"
%
% IK, to a distressed-looking 
% arts student at a dinner in
% Magdalen College

