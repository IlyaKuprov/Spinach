% Transforms the coefficients in front of Stevens operators, as
% produced by stevens.m, into the coefficients before the irredu-
% cible spherical tensor operators, as produced by irr_sph_ten.m
% function. Works up to 12th spherical rank. Source for ranks
% up to 6:
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
% Note: the scaling factors for ranks 7 to 12 were obtained by
%       projecting the operators returned by stevens.m onto the
%       operators returned by irr_sph_ten.m; the squared factors
%       are rational numbers that are the same for every spin
%       multiplicity and reproduce the published ranks 1 to 6.
%
% e.suturina@bath.ac.uk
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
a{7}=sqrt([32 16/7 416/7 104/7 4576/7 2288/7 13728/7 6864 13728/7 2288/7 4576/7 104/7 416/7 16/7 32]');
a{8}=sqrt([64 4 120 20/7 1040/7 156/7 1144/7 2860 823680 2860 1144/7 156/7 1040/7 20/7 120 4 64]');
a{9}=sqrt([128 64/9 2176/9 136/3 2720/9 272/63 7072/21 12771/65 570569/33 1555840 570569/33 12771/65 7072/21 272/63 2720/9 136/3 2176/9 64/9 128]');
a{10}=sqrt([256 64/5 2432/5 1216/15 82688/15 5168/3 10336/15 5168/15 537472/15 134368/5 11824384 134368/5 537472/15 5168/15 10336/15 5168/3 82688/15 1216/15 2432/5 64/5 256]');
a{11}=sqrt([512 256/11 10752/11 896/55 68096/55 34048/11 1157632/33 10336/33 82688/55 289408/55 2052166/3 22573824 2052166/3 289408/55 82688/55 10336/33 1157632/33 34048/11 68096/55 896/55 10752/11 256/11 512]');
a{12}=sqrt([1024 128/3 5888/3 2944/11 82432/33 20608/33 783104/11 55936/99 1613669/21 475456/11 950912/33 3328192/3 692263936 3328192/3 950912/33 475456/11 1613669/21 55936/99 783104/11 20608/33 82432/33 2944/11 5888/3 128/3 1024]');

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
   (~isscalar(k))||(mod(k,1)~=0)||(k<1)||(k>12)
    error('k must be a real integer between 1 and 12.');
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

