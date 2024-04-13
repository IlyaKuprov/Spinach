% Converts the nine components of the irreducible spherical tensor re-
% presentation of an interaction tensor into the Cartesian representa-
% tion with a 3x3 matrix. The conventions are matched to Equation (18)
% of the paper by Len Mueller (http://dx.doi.org/10.1002/cmr.a.20224).
%
% Spherical tensor components should be listed in the following order:
%
% rank 0: (0,0)
% rank 1: (1,1) (1,0) (1,-1)
% rank 2: (2,2) (2,1) (2,0) (2,-1) (2,-2)
%
% and should be supplied as coefficients in front of the corresponding 
% irreducible spherical tensor operators returned by irr_sph_ten.m fun-
% ction Syntax:
%
%                    M=sphten2mat(rank0,rank1,rank2)
%
% Parameters:
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
% Outputs:
%
%   M          - 3x3 interaction tensor
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=sphten2mat.m>

function M=sphten2mat(rank0,rank1,rank2)

% Check the input
grumble(rank0,rank1,rank2)

% Preallocate the answer
M=zeros(3);

% Rank 0 component
if ~isempty(rank0), M=M+rank0*eye(3); end

% Rank 1 components
if exist('rank1','var')&&~isempty(rank1)
    M=M-(1/2)*[0 0 -1; 0 0 -1i; 1 1i 0]*rank1(1);
    M=M+(1/sqrt(2))*[0 -1i 0; 1i 0 0; 0 0 0]*rank1(2);
    M=M-(1/2)*[0 0 -1; 0 0 1i; 1 -1i 0]*rank1(3);
end

% Rank 2 components
if exist('rank2','var')&&~isempty(rank2)
    M=M+(1/2)*[1 1i 0; 1i -1 0; 0 0 0]*rank2(1);
    M=M-(1/2)*[0 0 1; 0 0 1i; 1 1i 0]*rank2(2);
    M=M+(1/sqrt(6))*[-1 0 0; 0 -1 0; 0 0 2]*rank2(3);
    M=M+(1/2)*[0 0 1; 0 0 -1i; 1 -1i 0]*rank2(4);
    M=M+(1/2)*[1 -1i 0; -1i -1 0; 0 0 0]*rank2(5);
end

end

% Consistency enforcement
function grumble(rank0,rank1,rank2)
if (~isnumeric(rank0))||(~isnumeric(rank1))||(~isnumeric(rank2))
    error('all inputs must be vectors.');
end
if (numel(rank0)~=1)&&(~isempty(rank0))
    error('the first input must either be empty or have exactly one element.');
end
if (numel(rank1)~=3)&&(~isempty(rank1))
    error('the second input must either be empty or have exactly three elements.');
end
if (numel(rank2)~=5)&&(~isempty(rank2))
    error('the third input must either be empty or have exactly five elements.');
end
end

% When you're young, you look at television and think "there's a conspi-
% racy - the networks have conspired to dumb us down". But when you get
% a little older, you realize that's not true. The networks are in busi-
% ness to give people exactly what they want. That's a far more depres-
% sing thought.
%
% Steve Jobs

