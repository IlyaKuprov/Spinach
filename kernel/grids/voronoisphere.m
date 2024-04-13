% Voronoi tessellation on a sphere. Syntax:
%
%   [vertices,indices,polygons,sangles]=voronoisphere(xyz,res)
%
% Parameters:
%
%     xyz      - (3 x n) array, coordinates of n distinct vectors
%                 in R^3; these will be normalised
%
%     res      - polygon edges longer than this resolution will
%                have extra points added so that spherical plots
%                look neat (optional, default is pi/180)
%
% Outputs:
%
%     vertices - (3 x m) array, coordinates of the vertices of the
%                 Voronoi tessellation
%
%     indices  - (n x 1) cell array, j-th element contains the in-
%                 dices of the Voronoi cell vertices that corres-
%                 pond to xyz(:,j). Vertices are oriented counter-
%                 clockwise when looking from outside.
%
%     polygons - (n x 1) cell array, j-th element contains the dis-
%                 cretised spherical polygonal coordinates of the
%                 vertices of the j-th Voronoi cell.
%
%     sangles  - (n x 1) array, solid angles of each Voronoi cell
%
% brunoluong@yahoo.com
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=voronoisphere.m>

function [vertices,indices,polygons,sangles]=voronoisphere(xyz,res)

% Default polygon sampling resolution
if ~exist('res','var'), res=pi/180; end

% Check consistency and normalise
grumble(xyz); xyz=xyz./sqrt(sum(xyz.^2,1));

% Get Delaunay triangles of the convex hull
triangles=convhull(xyz.'); ntrian=size(triangles,1);

% Voronoi vertices are centres of Delaunay triangles
vertices=dtcentres(xyz,triangles);

% Compile triangle edge index
edges=sort([triangles(:,[1 2]); 
            triangles(:,[2 3]); 
            triangles(:,[3 1])],2);

% Link up forward and backward edge pairs
[~,~,J]=unique(edges,'rows');
if ~all(accumarray(J,1)==2)
    error('edge ambiguities found, check for point overlap.');
end

% Which 2 seeds the edge vertices correspond?
allids=repmat((1:ntrian).',[3 1]);
k=accumarray(J,(1:3*ntrian).',[],@(x){x});
k=[k{:}]; vid=allids(k.');

% each row is 2 cell ids of the edge
cellofedge = edges(k(1,:),:); % ne x 2
ne = size(cellofedge,1);
edges = repmat((1:ne).',[2 1]);
edgeofcell = accumarray(cellofedge(:),edges, [], @(x){x});

% Build the geodesic arcs that link two vertices
nverts = size(vid,1);
edgearcs = cell(1,nverts);
for k=1:nverts
    edgearcs{k} = geodark(vertices(:,vid(k,:)), res);
end

% Build the contour of the voronoi cells
polygons = cell(size(edgeofcell));
indices = cell(size(edgeofcell));
for k = 1:size(xyz,2)
    
    % Ordering and orientation of the edges
    v = cycling_edge(edgeofcell{k}, vid);
    v = oriented_edge(v, vertices, xyz(:,k));
    [~, loc] = ismember(sort(v,2),sort(vid,2),'rows');
    
    % Joining the arcs
    X = edgearcs(loc);
    flip = v(:,1)~=vid(loc,1);
    X(flip) = cellfun(@fliplr, X(flip), 'unif', false);
    X = cellfun(@(x) x(:,1:end-1), X, 'unif', false); % remove duplicated points
    polygons{k} = cat(2, X{:});
    
    % Indices of the Voronoi
    indices{k}=v(:,1);
    
end

% Return solid angles if needed
if nargout >= 4
    sangles = vcell_solidangle(vertices, indices, xyz);
end

end

% Returns a discretized arc between points A and B
function G = geodark(AB, resolution)
A = AB(:,1);
B = AB(:,2);
AxB = cross(A,B);
AdB = dot(A,B);
Ap = cross(AxB, A);
Ap = Ap/norm(Ap,2);
theta = atan2(sqrt(sum(AxB.^2,1)), AdB); % > 0
npnts = max(ceil(theta/resolution),2); % at least 2 points
theta = linspace(0, theta, npnts);
G = A*cos(theta) + Ap*sin(theta);
end

% Centres of Delaunay triangles
function cents=dtcentres(xyz,triangles)
xyz=reshape(xyz(:,triangles),[3 size(triangles)]);
A=xyz(:,:,1)-xyz(:,:,3); B=xyz(:,:,2)-xyz(:,:,3);
A2B=sum(A.^2,1).*B; B2A=sum(B.^2,1).*A; AxB=cross(A,B,1);
cents=xyz(:,:,3)+cross(A2B-B2A,AxB,1)./(2*sum(AxB.^2,1));
cents=cents.*sign(dot(AxB,xyz(:,:,3)))./sqrt(sum(cents.^2,1));
end

%%
function v = cycling_edge(edges, vertexes)
% Chain the edges in cycle
u = vertexes(edges,:).';
n = size(u, 2);
[~, ~, I] = unique(u);
I = reshape(I,[2 n]);
J = repmat(1:n,[2 1]);
if ~all(accumarray(I(:), 1) == 2)
    error('Topology issue due to numerical precision')
end
K = accumarray(I(:), J(:), [], @(x) {x});
K = [K{:}];
v = zeros([n 2]);
p = 0;
q = 1;
% chain the edges
for j = 1:n
    i = K(:,q);
    if i(1) == p
        p = i(2);
    else
        p = i(1);
    end    
    i = I(:,p);
    if i(1) == q
        v(j,:) = u([1 2],p);
        q = i(2);
    else
        v(j,:) = u([2 1],p);
        q = i(1);
    end
end % for-loop
end % cycling_edge

%%
function v = oriented_edge(v, P, xyz)
% Orient the edges counter-clockwise
Q = null(xyz.');
E = P(:,v([1:end 1],1));
xy = Q'*E;
a = (xy(1,1:end-1)-xy(1,2:end))*(xy(2,1:end-1)+xy(2,2:end))';
if xor(a < 0, det([xyz Q]) < 0) % Combine orientation and directness of Q
    v = rot90(v,2);
end
end % oriented_edge

% Consistency enforcement
function grumble(xyz)
if (~isnumeric(xyz))||(~isreal(xyz))||(size(xyz,1)~=3)
    error('xyz must be a [3 x N] array of real numbers.');
end
if size(xyz,2)<4, error('a minimum of four points needed.'); end
end

% May we never confuse honest dissent with disloyal subversion.
% 
% Dwight Eisenhower

