% Mesh preprocessing for drawing. Creates edge, triangle, and 
% rectangle data structures needed for fast plotting later.
% Syntax:
%                   mesh=mesh_preplot(mesh)
%
% Parameters:
%
%    mesh - Spinach mesh object
%
% Outputs:
%
%    mesh - updated mesh object
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=mesh_preplot.m>

function mesh=mesh_preplot(mesh)

% Check consistency
grumble(mesh);

% Prepare edge array for plotting 
nlines=size(mesh.idx.edges,1); A=zeros(1,3*nlines); B=zeros(1,3*nlines);
for n=1:nlines
    A((3*(n-1)+1):(3*n))=[mesh.x(mesh.idx.edges(n,1)) mesh.x(mesh.idx.edges(n,2)) NaN];
    B((3*(n-1)+1):(3*n))=[mesh.y(mesh.idx.edges(n,1)) mesh.y(mesh.idx.edges(n,2)) NaN];
end
mesh.plot.edg_a=A; mesh.plot.edg_b=B;

% Prepare triangle array for plotting
nlines=size(mesh.idx.triangles,1); A=zeros(1,5*nlines); B=zeros(1,5*nlines);
for n=1:nlines
    A((5*(n-1)+1):(5*n))=[mesh.x(mesh.idx.triangles(n,1)) mesh.x(mesh.idx.triangles(n,2)) ...
                          mesh.x(mesh.idx.triangles(n,3)) mesh.x(mesh.idx.triangles(n,1)) NaN];
    B((5*(n-1)+1):(5*n))=[mesh.y(mesh.idx.triangles(n,1)) mesh.y(mesh.idx.triangles(n,2)) ...
                          mesh.y(mesh.idx.triangles(n,3)) mesh.y(mesh.idx.triangles(n,1)) NaN];
end
mesh.plot.tri_a=A; mesh.plot.tri_b=B;

% Prepare rectangle array for plotting
nlines=size(mesh.idx.rectangles,1); A=zeros(1,6*nlines); B=zeros(1,6*nlines);
for n=1:nlines
    A((6*(n-1)+1):(6*n))=[mesh.x(mesh.idx.rectangles(n,1)) mesh.x(mesh.idx.rectangles(n,2)) ...
                          mesh.x(mesh.idx.rectangles(n,4)) mesh.x(mesh.idx.rectangles(n,3)) ...
                          mesh.x(mesh.idx.rectangles(n,1)) NaN];
    B((6*(n-1)+1):(6*n))=[mesh.y(mesh.idx.rectangles(n,1)) mesh.y(mesh.idx.rectangles(n,2)) ...
                          mesh.y(mesh.idx.rectangles(n,4)) mesh.y(mesh.idx.rectangles(n,3)) ...
                          mesh.y(mesh.idx.rectangles(n,1)) NaN];
end
mesh.plot.rec_a=A; mesh.plot.rec_b=B;

% Prepare Voronoi cell array for plotting
A=[]; B=[];
for n=1:numel(mesh.vor.cells)
    A=[A; mesh.vor.vertices(mesh.vor.cells{n},1); ...
          mesh.vor.vertices(mesh.vor.cells{n}(1),1); NaN]; %#ok<AGROW> 
    B=[B; mesh.vor.vertices(mesh.vor.cells{n},2); ...
          mesh.vor.vertices(mesh.vor.cells{n}(1),2); NaN]; %#ok<AGROW> 
end
mesh.plot.vor_a=A; mesh.plot.vor_b=B;

end

% Consistency enforcement
function grumble(mesh)
if ~isfield(mesh,'idx')
    error('vertex index is missing from spin_system.mesh structure.');
end
if ~isfield(mesh,'vor')
    error('Voronoi tessellation information is missing from spin_system.mesh structure.');
end
end

% Life is a tragedy for those who feel,
% and a comedy for those who think.
%
% Jean de la Bruyere

