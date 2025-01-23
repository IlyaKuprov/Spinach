% Plots a stick representation of a molecule from Cartesian coordinates
% supplied. Syntax:
%
%                          molplot(xyz,conmatrix)
%
% Parameters:
%
%             xyz - Cartesian coordinates, as Nx3 matrix, in
%                   Angstroms
%
%       conmatrix - NxN connectivity matrix indicating chemical
%                   bonds that should be drawn as sticks. If an
%                   empty vector is supplied, 1.6 Angstrom cut-
%                   off distance is used
%
% Outputs:
%
%       this function creates a figure
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=molplot.m>

function molplot(xyz,conmatrix)

% Check consistency
grumble(xyz,conmatrix);

% Get the connectivity matrix
if isempty(conmatrix)
    conmatrix=conmat(xyz,1.6);
end

% Prepare coordinate arrays
nbonds=nnz(conmatrix);
X=zeros(1,3*nbonds);
Y=zeros(1,3*nbonds);
Z=zeros(1,3*nbonds);
[rows,cols]=find(conmatrix);
for n=1:nbonds
    X((3*(n-1)+1):(3*n))=[xyz(rows(n),1) xyz(cols(n),1) NaN];
    Y((3*(n-1)+1):(3*n))=[xyz(rows(n),2) xyz(cols(n),2) NaN];
    Z((3*(n-1)+1):(3*n))=[xyz(rows(n),3) xyz(cols(n),3) NaN];
end

% Draw the molecule
plot3(X,Y,Z,'Color',[0.5 0.5 0.5],'LineWidth',1.5);

end

% Consistency enforcement
function grumble(xyz,conmatrix)
if (~isnumeric(xyz))||(~isreal(xyz))||(size(xyz,2)~=3)
    error('xyz must be an Nx3 real matrix.');
end
if ~isempty(conmatrix)
    if (size(conmatrix,1)~=size(conmatrix,2))||(size(conmatrix,1)~=size(xyz,1))
        error(['conmatrix must be a logical square matrix '...
               'of the same dimension as the number of rows in xyz.']);
    end
end
end

% The most terrifying fact about the universe is not that it is
% hostile but that it is indifferent... However vast the darkness,
% we must supply our own light.
% 
% Stanley Kubrick

