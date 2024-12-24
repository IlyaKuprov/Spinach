% Van Vleck perturbation theory, following Shavitt and Redmon, but
% excluding the quasi-degenerate split. Syntax:
%
%                     [Ep,G]=vvpert(E0,H1,order)
%
% Parameters:
%
%     E0    - eigenvalues of H0, a column vector of real 
%             numbers
%
%     H1    - perturbation, written in the basis that di-
%             agonalises H0
%
%     order - order of perturbation theory to be used, 5
%             is the maximum available
%
% Outputs:
%
%     Ep    - eigenvalues of H0+H1 to the specified order,
%             a column vector of reals, not necessarily 
%             sorted in the same way as the input
%
%     G     - Van Vleck transformation generator, such that
%             expm(G) is a square unitary matrix with eigen-
%             vectors in columns, in the same order as the 
%             eigenvalues in Ep
%
% Notes: there must be no degeneracies in H0; H1 must be Hermitian,
%        the theory only converges when norm(H1,2) is much than the
%        smallest energy gap in H0; complexity is linear in the or-
%        der and cubic in the matrix dimension.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=vvpert.m>

function [Ep,G]=vvpert(E0,H1,order)

% Check consistency
grumble(E0,H1,order);

% Reciprocal energy differences
Q=1./(E0'-E0); Q(logical(speye(size(Q))))=0;
if any(~isfinite(Q),'all')
    error('H0 has degenerate energy levels.');
end

% Split the perturbation
H1d=diag(diag(H1)); H1x=H1-H1d;

% Compute to specified order
G={}; W={};
if order>0
    G{1}=Q.*H1x;
    W{1}=H1d;
end
if order>1
    G{2}=Q.*   rocomm({H1d,G{1}});
    W{2}=(1/2)*rocomm({H1x,G{1}});
end
if order>2
    G{3}=Q.*(  rocomm({H1d,G{2}})+...
         (1/3)*rocomm({H1x,G{1},G{1}}));
    W{3}=(1/2)*rocomm({H1x,G{2}});
end
if order>3
    G{4}=Q.*(  rocomm({H1d,G{3}})+...
         (1/3)*rocomm({H1x,G{1},G{2}})+...
         (1/3)*rocomm({H1x,G{2},G{1}}));
    W{4}=(1/2)*rocomm({H1x,G{3}})-...
        (1/24)*rocomm({H1x,G{1},G{1},G{1}});
end
if order>4
    G{5}=Q.*(  rocomm({H1d,G{4}})+...
         (1/3)*rocomm({H1x,G{1},G{3}})+...
         (1/3)*rocomm({H1x,G{2},G{2}})+...
         (1/3)*rocomm({H1x,G{3},G{1}})-...
        (1/45)*rocomm({H1x,G{1},G{1},G{1},G{1}}));
    W{5}=(1/2)*rocomm({H1x,G{4}})-...
        (1/24)*rocomm({H1x,G{1},G{1},G{2}})-...
        (1/24)*rocomm({H1x,G{1},G{2},G{1}})-...
        (1/24)*rocomm({H1x,G{2},G{1},G{1}});
end
if order>5
    error('VVPT not implemented for orders higher than five.');
end

% Get the energies
W=cellfun(@diag,W,'UniformOutput',0);
Ep=E0+sum(cell2mat(W),2);

% Get the generator
G=reshape(G,[1 1 numel(G)]);
G=sum(cell2mat(G),3);

end

% Consistency enforcement
function grumble(E0,H1,order)
if (~isnumeric(E0))||(~isreal(E0))||(~iscolumn(E0))
    error('E0 must be a real column vector.');
end
if (~isnumeric(H1))||(size(H1,1)~=numel(E0))||...
   (size(H1,2)~=numel(E0))
    error('H1 must be a Hermitian matrix with dimensions matching E0');
end
if (~isnumeric(order))||(~isreal(order))||(~isscalar(order))||...
   (mod(order,1)~=0)||(order<1)||(order>5)
    error('order must be a real integer between 1 and 5');
end
end

% When have we ever had stability and security in this
% region that we should be concerned about losing it?
%
% Kurdish officials, to BBC

