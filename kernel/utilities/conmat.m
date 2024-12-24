% Molecular connectivity matrix calculator with N*log(N) 
% asymptotic complexity scaling with respect to the num-
% ber or atoms. Syntax:
%
%                conmatrix=conmat(xyz,r0)
%
% Parameters:
%
%       xyz  - an array with N rows and three columns,
%              giving the Cartesian coordinates of
%              each particle
%
%        r0  - the distance below which the particles 
%              are to be considered "connected"
%
% Output:
%
%  conmatrix - a sparse logical matrix containing 1
%              at the positions corresponding to the 
%              connected particles
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=conmat.m>

function conmatrix=conmat(xyz,r0)

% Check consistency
grumble(xyz,r0);

% Sort X coordinates
[x_sorted,x_index]=sort(xyz(:,1));

% Scan X coordinates
A=false(size(xyz,1),size(xyz,1));
for n=1:size(xyz,1)
    for k=(n+1):size(xyz,1)
        if abs(x_sorted(n)-x_sorted(k))<r0
            A(x_index(n),x_index(k))=1; 
            A(x_index(k),x_index(n))=1;
        else
            break;
        end
    end
end

% Sort Y coordinates
[y_sorted,y_index]=sort(xyz(:,2));

% Scan Y coordinates
B=false(size(xyz,1),size(xyz,1));
for n=1:size(xyz,1)
    for k=(n+1):size(xyz,1)
        if abs(y_sorted(n)-y_sorted(k))<r0
            B(y_index(n),y_index(k))=1; 
            B(y_index(k),y_index(n))=1;
        else
            break;
        end
    end
end

% Sort Z coordinates
[z_sorted,z_index]=sort(xyz(:,3));

% Scan Y coordinates
C=false(size(xyz,1),size(xyz,1));
for n=1:size(xyz,1)
    for k=(n+1):size(xyz,1)
        if abs(z_sorted(n)-z_sorted(k))<r0
            C(z_index(n),z_index(k))=1; 
            C(z_index(k),z_index(n))=1;
        else
            break;
        end
    end
end

% Compile the dirty matrix
conmatrix=A&B&C;

% Clean up the matrix
[row,col]=find(conmatrix);
for n=1:numel(row)
    if norm(xyz(row(n),:)-xyz(col(n),:),2)>r0
        conmatrix(row(n),col(n))=0;
    end
end
conmatrix=sparse(conmatrix);
        
end

% Consistency enforcement
function grumble(xyz,r0)
if (~isnumeric(xyz))||(size(xyz,2)~=3)||(~isreal(xyz))
    error('xyz parameter should be a real matrix with three columns.');
end
if (~isnumeric(r0))||(~isreal(r0))||(r0<=0)
    error('r0 must be a positive real number.');
end
end

% It's always good to be underestimated.
%
% Donald Trump

