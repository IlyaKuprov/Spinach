% Testing IK's favourite equation numerically - auxiliary
% matrix expression against a high-accuracy finite diffe-
% rence approximation.
%
% i.kuprov@soton.ac.uk

function auxmat_test()

% Random matrices, arbitrary function
dim=randi(20); f=@(A)expm(A);
A=randn(dim)+1i*randn(dim); A=A/norm(A,2);
B=randn(dim)+1i*randn(dim); B=B/norm(B,2);
C=randn(dim)+1i*randn(dim); C=C/norm(B,2);

% Some weird ass function and its derivatives
L=@(a,b,c)(a*A+c*cos(b)*B+(a*c^4)*C+a*b*c*B*C);
dL_da=@(a,b,c)(A+(c^4)*C+b*c*B*C);
dL_db=@(a,b,c)(-c*sin(b)*B+a*c*B*C);
dL_dc=@(a,b,c)(cos(b)*B+(4*a*c^3)*C+a*b*B*C);

% Real coefficients and fin. diff. increment
x=randn(1); y=randn(1); z=randn(1); h=sqrt(eps);

% First parameter
P=f([  L(x,y,z)   dL_da(x,y,z);
     zeros(dim)       L(x,y,z)]);
P=P(1:dim,(dim+1):(2*dim));
Q=2*(f(L(x+h,y,z))-f(L(x-h,y,z)))/(3*h)-...
    (f(L(x+2*h,y,z))-f(L(x-2*h,y,z)))/(12*h);
if norm(P-Q,2)/norm(P,2)>1e-3
    error('differentiation test failed');
else
    disp('differentiation test passed');
end

% Second parameter
P=f([  L(x,y,z)   dL_db(x,y,z);
     zeros(dim)       L(x,y,z)]);
P=P(1:dim,(dim+1):(2*dim));
Q=2*(f(L(x,y+h,z))-f(L(x,y-h,z)))/(3*h)-...
    (f(L(x,y+2*h,z))-f(L(x,y-2*h,z)))/(12*h);
if norm(P-Q,2)/norm(P,2)>1e-3
    error('differentiation test failed');
else
    disp('differentiation test passed');
end

% Third parameter
P=f([  L(x,y,z)   dL_dc(x,y,z);
     zeros(dim)       L(x,y,z)]);
P=P(1:dim,(dim+1):(2*dim));
Q=2*(f(L(x,y,z+h))-f(L(x,y,z-h)))/(3*h)-...
    (f(L(x,y,z+2*h))-f(L(x,y,z-2*h)))/(12*h);
if norm(P-Q,2)/norm(P,2)>1e-3
    error('differentiation test failed');
else
    disp('differentiation test passed');
end

end

