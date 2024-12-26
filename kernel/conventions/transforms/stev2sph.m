% Transforms the coefficients in front of Stevens operators, as 
% produced by stevens.m, into the coefficients before the irredu-
% cible spherical tensor operators, as produced by irr_sph_ten.m
% function. Works up to 6th spherical rank. Source:
%
%         http://dx.doi.org/10.1088/0022-3719/18/7/009
% 
% Syntax:
%
%                       Bkq=stev2sph(k,Bkq)
%
% Parameters:
%
%   k     - the spherical rank in question
%
%   Bkq   - a column of 2k+1 real coefficients
%           in front of Stevens operators, in
%           increasing order of projections
%
% Outputs:
%
%   Bkq   - a column of 2k+1 complex coefficients
%           in front of irreducible spherical 
%           tensor operators, in decreasing order
%           of projections
%
% e.suturina@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=stev2sph.m>

function Bkq=stev2sph(k,Bkq)

% Check consistency
grumble(k,Bkq);

% Catalog the stupid scaling factors
a{1}=[1/sqrt(2) 1 1/sqrt(2)]';
a{2}=[1 1/2 sqrt(6) 1/2 1]';
a{3}=[sqrt(2) 1/sqrt(3) sqrt(10/3) sqrt(10) sqrt(10/3) 1/sqrt(3) sqrt(2)]';
a{4}=[2 1/sqrt(2) sqrt(7) sqrt(7/2) 2*sqrt(70) sqrt(7/2) sqrt(7) 1/sqrt(2) 2]';
a{5}=[2*sqrt(2) 2/sqrt(5) 6*sqrt(2/5) sqrt(3/5) 2*sqrt(21/5) 6*sqrt(14) 2*sqrt(21/5) sqrt(3/5) 6*sqrt(2/5) 2/sqrt(5) 2*sqrt(2)]';
a{6}=[4 2/sqrt(3) 4*sqrt(11/6) 2*sqrt(11/5) 4*sqrt(11/5) sqrt(22) 4*sqrt(231) sqrt(22) 4*sqrt(11/5) 2*sqrt(11/5) 4*sqrt(11/6) 2/sqrt(3) 4]';

% Form the transformation matrix diagonal
criss=[-1i*(-1).^(k:-1:1)'; 1; ones(k,1)   ].*a{k};

% Form the transformation matrix antidiagonal
cross=[+1i*ones(k,1);       0; (-1).^(1:k)'].*a{k};        

% Form the transformation matrix 
A=diag(criss)+fliplr(diag(cross));

% Transform the coefficients
Bkq=transpose(Bkq'*A);

end

% Consistency enforcement
function grumble(k,Bkq)
if (~isnumeric(k))||(~isreal(k))||(~isfinite(k))||...
   (~isscalar(k))||(mod(k,1)~=0)||(k<1)||(k>6)
    error('k must be a real integer between 1 and 6.');
end
if (~isnumeric(Bkq))||(~isreal(Bkq))||any(~isfinite(Bkq))||...
   (~iscolumn(Bkq))||(numel(Bkq)~=2*k+1)
    error('Bkq must be a column vector with 2*k+1 real elements.');
end
end

% K.W.H. Stevens has done a great disservice to Magnetic Resonance by
% his ill-considered choice of basis operators for the crystal field 
% theory he was developing. Choosing irreducible spherical tensors in-
% stead would have saved many days to nearly everybody in this field.
% At the moment, the community is stuck with the ridiculously bad de-
% finitions that do not follow the 3D rotation group - "for historical
% reasons". Someone would have created crystal field theory if Stevens
% hadn't - but he had poisoned it forever by doing it badly.

