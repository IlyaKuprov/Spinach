% 2D microfluidic mesh cropping. Updates the mesh object to 
% remove anything outside the user-specified vertex coordi-
% nate ranges. Syntax:
%
%     spin_system=mesh_crop(spin_system,xrange,yrange)
%
% Parameters:
%
%    spin_system - Spinach data structure with 
%                  a .mesh subfield present
%
%    xrange      - two-element vector specify-
%                  ing X axis range, [min max]
%
%    yrange      - two-element vector specify-
%                  ing Y axis range, [min max]
%
% Outputs:
%
%    spin_system - updated data structure
%
% a.acharya@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=mesh_crop.m>

function spin_system=mesh_crop(spin_system,xrange,yrange)

% Check consistency
grumble(spin_system,xrange,yrange);

% Remove tessellation and preplot
if isfield(spin_system.mesh,'vor')
    spin_system.mesh=rmfield(spin_system.mesh,'vor');
    report(spin_system,'WARNING: out-of-date tessellation information removed.');
end
if isfield(spin_system.mesh,'plot')
    spin_system.mesh=rmfield(spin_system.mesh,'plot');
    report(spin_system,'WARNING: out-of-date preplot information removed.');
end

% Find vertices in the user-specified range
vertices_in_range=(spin_system.mesh.x>=xrange(1))&(spin_system.mesh.x<=xrange(2))&...
                  (spin_system.mesh.y>=yrange(1))&(spin_system.mesh.y<=yrange(2));
vertices_in_range=find(vertices_in_range);

% Find edges in the user-specified range
edges_in_range=ismember(spin_system.mesh.idx.edges(:,1),vertices_in_range)&...
               ismember(spin_system.mesh.idx.edges(:,2),vertices_in_range);
spin_system.mesh.idx.edges=spin_system.mesh.idx.edges(edges_in_range,:);

% Re-index edges with updated vertices 
[~,spin_system.mesh.idx.edges(:,1)]=ismember(spin_system.mesh.idx.edges(:,1),vertices_in_range); 
[~,spin_system.mesh.idx.edges(:,2)]=ismember(spin_system.mesh.idx.edges(:,2),vertices_in_range);

% Find triangles in the user-specified range 
tris_in_range=ismember(spin_system.mesh.idx.triangles(:,1),vertices_in_range)&...
              ismember(spin_system.mesh.idx.triangles(:,2),vertices_in_range)&...
              ismember(spin_system.mesh.idx.triangles(:,3),vertices_in_range);
spin_system.mesh.idx.triangles=spin_system.mesh.idx.triangles(tris_in_range,:);

% Re-index trianges with updated vertices 
[~,spin_system.mesh.idx.triangles(:,1)]=ismember(spin_system.mesh.idx.triangles(:,1),vertices_in_range);
[~,spin_system.mesh.idx.triangles(:,2)]=ismember(spin_system.mesh.idx.triangles(:,2),vertices_in_range);
[~,spin_system.mesh.idx.triangles(:,3)]=ismember(spin_system.mesh.idx.triangles(:,3),vertices_in_range);

% Find rectangles in the user-specified range
rect_in_range=ismember(spin_system.mesh.idx.rectangles(:,1),vertices_in_range)&...
              ismember(spin_system.mesh.idx.rectangles(:,2),vertices_in_range)&...
              ismember(spin_system.mesh.idx.rectangles(:,3),vertices_in_range)&...
              ismember(spin_system.mesh.idx.rectangles(:,4),vertices_in_range);
spin_system.mesh.idx.rectangles=spin_system.mesh.idx.rectangles(rect_in_range,:);

% Re-index rectangles with updated vertices 
[~,spin_system.mesh.idx.rectangles(:,1)]=ismember(spin_system.mesh.idx.rectangles(:,1),vertices_in_range);     
[~,spin_system.mesh.idx.rectangles(:,2)]=ismember(spin_system.mesh.idx.rectangles(:,2),vertices_in_range); 
[~,spin_system.mesh.idx.rectangles(:,3)]=ismember(spin_system.mesh.idx.rectangles(:,3),vertices_in_range);
[~,spin_system.mesh.idx.rectangles(:,4)]=ismember(spin_system.mesh.idx.rectangles(:,4),vertices_in_range);

% Crop coordinates, velocities and concentrations
spin_system.mesh.x=spin_system.mesh.x(vertices_in_range); 
spin_system.mesh.y=spin_system.mesh.y(vertices_in_range);  
if isfield(spin_system.mesh,'u'), spin_system.mesh.u=spin_system.mesh.u(vertices_in_range); end
if isfield(spin_system.mesh,'v'), spin_system.mesh.v=spin_system.mesh.v(vertices_in_range); end
if isfield(spin_system.mesh,'c'), spin_system.mesh.c=spin_system.mesh.c(vertices_in_range,:); end

% First guess active vertices
if isfield(spin_system.mesh.idx,'active')
    report(spin_system,'WARNING: active vertex list overwritten.');
end
spin_system.mesh.idx.active=unique(spin_system.mesh.idx.triangles(:));

end

% Consistency enforcement
function grumble(spin_system,xrange,yrange)
if ~isfield(spin_system,'mesh')
    error('mesh information is missing from the spin_system structure.');
end
if ~isfield(spin_system.mesh,'idx')
    error('vertex index is missing from spin_system.mesh structure.');
end
if (~isnumeric(xrange))||(~isreal(xrange))||(numel(xrange)~=2)||(xrange(1)>=xrange(2))
    error('xrange must be a two-element real vector in ascending order.');
end
if (~isnumeric(yrange))||(~isreal(yrange))||(numel(yrange)~=2)||(yrange(1)>=yrange(2))
    error('yrange must be a two-element real vector in ascending order.');
end
end

% I will not utter falsehoods, but have no objection to 
% meaningless statements.
%
% A.J. Ayer, a noted atheist, about saying
% a grace at the New College High Table

