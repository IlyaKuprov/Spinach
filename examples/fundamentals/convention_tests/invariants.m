% A test of Equation 3 in http://dx.doi.org/10.1002/chem.200902300
%
% i.kuprov@soton.ac.uk

function invariants()

% Rhombic test
for n=1:100
    
    % Random symmetric matrix
    A=randn(3); A=A+A';
    
    % Eigenvalues
    [~,D]=eig(A); D=diag(D); D=D(randperm(3));
    
    % Axiality and rhombicity
    Ax=2*D(3)-(D(1)+D(2)); Rh=D(1)-D(2);
    
    % Standard relaxation theory invariant
    P=(Ax^2+3*Rh^2)/6;
    
    % Neat eigenvalue form
    Q=(2/3)*(D(1)^2+D(2)^2+D(3)^2-D(1)*D(2)-D(1)*D(3)-D(2)*D(3));
    
    % Difference and diagnostics
    if abs(P-Q)>1e-6, error('Rhombic test failed.'); end
    
end

% Successful completion message
if n==100, disp('Rhombic test passed.'); end

% Axial test
for n=1:100
    
    % Random axial eigenvalues
    D=randn(2,1); D(3)=D(2); 
    
    % Zero trace
    D=D-sum(D)/3;
    
    % Standard invarinat
    P=(D(1)-D(2))^2;
    
    % Random shuffle
    D=D(randperm(3));
    
    % Neat invariant
    Q=D(1)^2+D(2)^2+D(3)^2-D(1)*D(2)-D(1)*D(3)-D(2)*D(3);
    
    % Difference and diagnostics
    if abs(P-Q)>1e-6, error('Axial test failed.'); end
    
end

% Successful completion message
if n==100, disp('Axial test passed.'); end

end

