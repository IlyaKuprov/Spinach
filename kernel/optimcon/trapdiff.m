% Directional derivatives for the trapezium product quadrature publi-
% shed by Iserles and Norsett (see Corollary 3.3) in
%
%               https://doi.org/10.1098/rsta.1999.0362
% 
% The derivatives are of the following propagator:
%
%          expm(-i*((HL+HR)/2+i*dt*(sqrt(3)/12)*[HL,HR])*dt)
% 
% with respect to the coefficients cL,cR in the evolution generators
% HL and HR on the left and the right side of the interval respecti-
% vely. Evolution generators HL and HR are split into the drift part
% Ho and the control part Hc, such that HL=Ho+cL*Hc and HR=Ho+cR*Hc
% on the left and the right edge of the interval. 
%
% The derivatives are calculated using Eq 16 of Goodwin and Kuprov
%
%                 https://doi.org/10.1063/1.4928978
%
% Syntax:
%
%             [DL,DR]=trapdiff(spin_system,Hd,Hc,dt,cL,cR)
%
% Parameters:
%
%     Hd - a cell array of two matrices containing drift 
%          generators at the left (first element) and the
%          right (second element) edge of the interval
%
%     Hc - control operator or superoperator
%         
%     dt - interval duration, seconds
%
%     cL - control operator coefficient at the
%          left edge of the interval
%
%     cR - control operator coefficient at the
%          right edge of the interval
%
% Outputs:
%
%     DL - derivative of the interval propagator
%          with respect to cL
%
%     DR - derivative of the interval propagator
%          with respect to cR
%
% u.rasulov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=trapdiff.m>

function [DL,DR]=trapdiff(spin_system,Hd,Hc,dt,cL,cR)

% Check consistency
grumble(Hd,Hc,dt,cL,cR)

% Precompute directions
H_dir_L=(1/2)*Hc+1i*dt*(sqrt(3)/12)*(Hc*Hd{2}-Hd{2}*Hc);
H_dir_R=(1/2)*Hc+1i*dt*(sqrt(3)/12)*(Hd{1}*Hc-Hc*Hd{1});

% Call directional derivative function
D=dirdiff(spin_system,(Hd{1}+Hd{2})/2+cL*H_dir_L+cR*H_dir_R,H_dir_L,dt,2); DL=D{2};
D=dirdiff(spin_system,(Hd{1}+Hd{2})/2+cL*H_dir_L+cR*H_dir_R,H_dir_R,dt,2); DR=D{2};

end

% Consistency enforcement
function grumble(Hd,Hc,dt,cL,cR)
if ~iscell(Hd), error('Hd must be a cell array of matrices.'); end
for n=1:numel(Hd)     
    if (~isnumeric(Hd{n}))||(size(Hd{n},1)~=size(Hd{n},2))
        error('all elements of Hd cell array must be square matrices.');
    end
end
if (~isnumeric(Hc))||(size(Hc,1)~=size(Hc,2))||(size(Hc,1)~=size(Hd{1},1))
    error('Hc must have the same size as drift generators.');
end

if (~isnumeric(dt))||(~isreal(dt))||(~isscalar(dt))||(dt<=0)
    error('dt must be a positive real scalar.');
end
if (~isnumeric(cL))||(~isreal(cL))||(~isnumeric(cR))||(~isreal(cR))
    error('cL and cR must be real scalars.');
end
end

% Don't get mad at the mirror when
% your mug's askew. 
% 
% Vladimir Putin

