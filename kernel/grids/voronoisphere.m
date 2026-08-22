% Voronoi tessellation of the unit sphere around the specified po-
% ints, computed as the exact geometric dual of the Delaunay tri-
% angulation returned by the convex hull. Syntax:
%
%      [vertices,indices,polygons,sangles]=voronoisphere(xyz)
%
% Parameters:
%
%     xyz      - (3 x n) array, coordinates of n distinct vectors
%                in R^3; these will be normalised
%
% Outputs:
%
%     vertices - (3 x m) array, coordinates of the vertices of the
%                Voronoi tessellation
%
%     indices  - (n x 1) cell array, j-th element contains the in-
%                dices of the Voronoi cell vertices that correspond
%                to xyz(:,j). Vertices are oriented counterclockwise
%                when looking from outside.
%
%     polygons - (n x 1) cell array, j-th element contains the coor-
%                dinates of the vertices of the j-th Voronoi cell,
%                in the same counterclockwise order
%
%     sangles  - (n x 1) array, solid angles of each Voronoi cell
%
% Note: every Voronoi vertex is the circumcentre of a Delaunay tri-
%       angle, placed on the side of the triangle plane that the
%       convex hull property guarantees to be empty; every Voronoi
%       cell is geodesically convex and therefore star-shaped around
%       its generator point, which makes the angular sort used here
%       an exact construction rather than a heuristic.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=voronoisphere.m>

function [vertices,indices,polygons,sangles]=voronoisphere(xyz)

% Check consistency
grumble(xyz);

% Normalise the points
xyz=xyz./sqrt(sum(xyz.^2,1)); npts=size(xyz,2);

% Delaunay triangles of a spherical point set are convex hull facets
triangles=convhull(xyz.'); ntrian=size(triangles,1);

% Refuse point sets whose hull is not a closed surface
edges=sort([triangles(:,[1 2]); triangles(:,[2 3]); triangles(:,[3 1])],2);
[~,~,edge_ids]=unique(edges,'rows');
if ~all(accumarray(edge_ids,1)==2)
    error('the point set does not triangulate into a closed surface.');
end

% Voronoi vertices are triangle circumcentres on the empty side
normals=cross(xyz(:,triangles(:,2))-xyz(:,triangles(:,1)),...
              xyz(:,triangles(:,3))-xyz(:,triangles(:,1)),1);
vertices=normals./sqrt(sum(normals.^2,1));
empty_side=sign(sum(vertices.*(xyz(:,triangles(:,1))-mean(xyz,2)),1));
vertices=vertices.*empty_side;

% List the triangles incident at each point
incidence=accumarray([triangles(:,1); triangles(:,2); triangles(:,3)],...
                     repmat((1:ntrian).',3,1),[npts 1],@(x){x});

% Build each cell by angular sorting around its generator
indices=cell(npts,1); polygons=cell(npts,1);
for n=1:npts

    % Right-handed tangent frame at the generator
    pivot=zeros(3,1); [~,pos]=min(abs(xyz(:,n))); pivot(pos)=1;
    tang_a=cross(xyz(:,n),pivot); tang_a=tang_a/norm(tang_a,2);
    tang_b=cross(xyz(:,n),tang_a);

    % Sort the cell vertices counterclockwise from outside
    cell_verts=vertices(:,incidence{n});
    [~,ccw_order]=sort(atan2(tang_b.'*cell_verts,tang_a.'*cell_verts));
    indices{n}=incidence{n}(ccw_order);
    polygons{n}=vertices(:,indices{n});

end

% Return solid angles if needed
if nargout>=4
    sangles=vcell_solidangle(vertices,indices,xyz);
end

end

% Consistency enforcement
function grumble(xyz)
if (~isnumeric(xyz))||(~isreal(xyz))||(size(xyz,1)~=3)
    error('xyz must be a [3 x N] array of real numbers.');
end
if any(~isfinite(xyz(:)))
    error('xyz must not contain NaN or Inf elements.');
end
if size(xyz,2)<4
    error('a minimum of four points is needed.');
end
if any(sum(xyz.^2,1)==0)
    error('xyz must not contain zero-length vectors.');
end
xyz=xyz./sqrt(sum(xyz.^2,1));
if size(uniquetol(xyz.',1e-9,'ByRows',true),1)<size(xyz,2)
    error('coincident points detected on the unit sphere.');
end
if min(svd(xyz))<1e-12*size(xyz,2)
    error('the points must not all lie on a single great circle.');
end
end

% May we never confuse honest dissent with disloyal subversion.
%
% Dwight Eisenhower

