% Test of the span-skew interaction convention.
%
% i.kuprov@soton.ac.uk

function spsk_test()

% Eigenvalues
xx=rand(); 
yy=rand()+3; 
zz=rand()+6;

% Euler angles
alp=2*pi*rand(); 
bet=pi*rand(); 
gam=2*pi*rand();

% Manual construction
R=euler2dcm(alp,bet,gam);
AM=R*diag([xx yy zz])*R';

% Spinach construction
iso=(xx+yy+zz)/3; 
sp=zz-xx; sk=3*(yy-iso)/sp;
AS=spsk2mat(iso,sp,sk,alp,bet,gam);

% Difference
if norm(AM-AS,1)<1e-6
    disp('Test passed.');
else
    error('Test FAILED.');
end

end

