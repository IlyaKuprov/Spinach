% Converts a 3x3 interaction matrix into the irreducible spherical tensor
% notation: one rank 0 component, three rank 1 components and five rank 2
% components to the total of nine independent components. The conventions
% are matched to Equation (22) of the paper by Len Mueller:
%
%                 http://dx.doi.org/10.1002/cmr.a.20224
%
% The components are listed in the following order:
%
% rank 0: (0,0)
% rank 1: (1,1) (1,0) (1,-1)
% rank 2: (2,2) (2,1) (2,0) (2,-1) (2,-2)
%
% and are returned as coefficients in front of the corresponding irredu-
% cible spherical tensor operators returned by irr_sph_ten.m function.
% Syntax:
%
%                    [rank0,rank1,rank2]=mat2sphten(M)
%
% Parameters:
%
%   M          - 3x3 interaction tensor
%
% Outputs:
%
%   rank0      - a single number giving the coefficient of T(0,0) in
%                the spherical tensor expansion.
%
%   rank1      - a row vector with three numbers giving the coeffici-
%                ents of T(1,1), T(1,0) and T(1,-1) in the spherical
%                tensor expansion.
%
%   rank2      - a row vector with five numbers giving the coeffici-
%                ents of T(2,2), T(2,1), T(2,0), T(2,-1) and T(2,-2)
%                in the spherical tensor expansion.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=mat2sphten.m>

function [rank0,rank1,rank2]=mat2sphten(M)

% Adapt to empty matrices
if isempty(M), M=zeros(3); end

% Check consistency
grumble(M);

% Set result dimensions
rank1=zeros(3,1); rank2=zeros(5,1);

% Rank 0 component
rank0=trace(M)/3;

% Rank 1 components
rank1(1)=-(1/2)*(M(3,1)-M(1,3)-1i*(M(3,2)-M(2,3)));
rank1(2)=+(1i/sqrt(2))*(M(1,2)-M(2,1));
rank1(3)=-(1/2)*(M(3,1)-M(1,3)+1i*(M(3,2)-M(2,3)));

% Rank 2 components
rank2(1)=+(1/2)*(M(1,1)-M(2,2)-1i*(M(1,2)+M(2,1)));
rank2(2)=-(1/2)*(M(1,3)+M(3,1)-1i*(M(2,3)+M(3,2)));
rank2(3)=+(1/sqrt(6))*(2*M(3,3)-M(1,1)-M(2,2));
rank2(4)=+(1/2)*(M(1,3)+M(3,1)+1i*(M(2,3)+M(3,2)));
rank2(5)=+(1/2)*(M(1,1)-M(2,2)+1i*(M(1,2)+M(2,1)));

end

% Consistency enforcement
function grumble(M)
if (~isnumeric(M))||(~isreal(M))||...
   (~ismatrix(M))||any(size(M)~=[3 3])
    error('the argument must be a real 3x3 matrix.');
end
end

% When I was 13 I think - Hewlett and Packard were my idols - I called up
% Bill Hewlett because he lived in Palo Alto and there were no unlisted
% numbers in the phonebook - which gives you a clue to my age. And he
% picked up the phone and I talked to him and I asked him if he'd give me
% some spare parts for something I was building called a frequency coun-
% ter. And he did, but in addition to that he gave me something way more
% important. He gave me a job that summer - a summer job - at Hewlett-
% Packard right here in Santa Clara off 280, in a division that built fre-
% quency counters. And I was in heaven.
%
% Steve Jobs

