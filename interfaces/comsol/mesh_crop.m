% 2D microfluidic mesh cropping. Updates the mesh object to 
% remove anything outside the user-specified vertex coordi-
% nate ranges. Syntax:
%
%               mesh=mesh_crop(mesh,ranges)
%
% Parameters:
%
%    mesh        - Spinach mesh object
%
%    ranges      - {[xmin xmax],[ymin ymax]}
%
% Outputs:
%
%    mesh        - updated mesh object
%
% a.acharya@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=mesh_crop.m>

function mesh=mesh_crop(mesh,ranges)

% Check consistency
grumble(mesh,ranges);

% Remove tessellation and preplot
if isfield(mesh,'vor')
    mesh=rmfield(mesh,'vor');
    report(spin_system,'WARNING: out-of-date tessellation information removed.');
end
if isfield(mesh,'plot')
    mesh=rmfield(mesh,'plot');
    report(spin_system,'WARNING: out-of-date preplot information removed.');
end

% Find vertices in the user-specified range
vertices_in_range=(mesh.x>=ranges{1}(1))&(mesh.x<=ranges{1}(2))&...
                  (mesh.y>=ranges{2}(1))&(mesh.y<=ranges{2}(2));
vertices_in_range=find(vertices_in_range);

% Find edges in the user-specified range
edges_in_range=ismember(mesh.idx.edges(:,1),vertices_in_range)&...
               ismember(mesh.idx.edges(:,2),vertices_in_range);
mesh.idx.edges=mesh.idx.edges(edges_in_range,:);

% Re-index edges with updated vertices 
[~,mesh.idx.edges(:,1)]=ismember(mesh.idx.edges(:,1),vertices_in_range); 
[~,mesh.idx.edges(:,2)]=ismember(mesh.idx.edges(:,2),vertices_in_range);

% Find triangles in the user-specified range 
tris_in_range=ismember(mesh.idx.triangles(:,1),vertices_in_range)&...
              ismember(mesh.idx.triangles(:,2),vertices_in_range)&...
              ismember(mesh.idx.triangles(:,3),vertices_in_range);
mesh.idx.triangles=mesh.idx.triangles(tris_in_range,:);

% Re-index trianges with updated vertices 
[~,mesh.idx.triangles(:,1)]=ismember(mesh.idx.triangles(:,1),vertices_in_range);
[~,mesh.idx.triangles(:,2)]=ismember(mesh.idx.triangles(:,2),vertices_in_range);
[~,mesh.idx.triangles(:,3)]=ismember(mesh.idx.triangles(:,3),vertices_in_range);

% Find rectangles in the user-specified range
rect_in_range=ismember(mesh.idx.rectangles(:,1),vertices_in_range)&...
              ismember(mesh.idx.rectangles(:,2),vertices_in_range)&...
              ismember(mesh.idx.rectangles(:,3),vertices_in_range)&...
              ismember(mesh.idx.rectangles(:,4),vertices_in_range);
mesh.idx.rectangles=mesh.idx.rectangles(rect_in_range,:);

% Re-index rectangles with updated vertices 
[~,mesh.idx.rectangles(:,1)]=ismember(mesh.idx.rectangles(:,1),vertices_in_range);     
[~,mesh.idx.rectangles(:,2)]=ismember(mesh.idx.rectangles(:,2),vertices_in_range); 
[~,mesh.idx.rectangles(:,3)]=ismember(mesh.idx.rectangles(:,3),vertices_in_range);
[~,mesh.idx.rectangles(:,4)]=ismember(mesh.idx.rectangles(:,4),vertices_in_range);

% Crop coordinates, velocities and concentrations
mesh.x=mesh.x(vertices_in_range); 
mesh.y=mesh.y(vertices_in_range);  
if isfield(mesh,'u'), mesh.u=mesh.u(vertices_in_range); end
if isfield(mesh,'v'), mesh.v=mesh.v(vertices_in_range); end
if isfield(mesh,'c'), mesh.c=mesh.c(vertices_in_range,:); end

% First guess active vertices
if isfield(mesh.idx,'active')
    report(spin_system,'WARNING: active vertex list was overwritten.');
end
mesh.idx.active=unique(mesh.idx.triangles(:));

end

% Consistency enforcement
function grumble(mesh,ranges)
if ~isfield(mesh,'idx')
    error('vertex index is missing from the mesh structure.');
end
if (~iscell(ranges))||(numel(ranges)~=2)||...
   (~isnumeric(ranges{1}))||(~isreal(ranges{1}))||...
   (numel(ranges{1})~=2)||(ranges{1}(1)>=ranges{1}(2))||...
   (~isnumeric(ranges{2}))||(~isreal(ranges{2}))||...
   (numel(ranges{2})~=2)||(ranges{2}(1)>=ranges{2}(2))
    error('ranges must be {[xmin xmax],[ymin ymax]}');
end
end

% I will not utter falsehoods, but have no objection to 
% meaningless statements.
%
% A.J. Ayer, a noted atheist, about saying
% a grace at the New College High Table

