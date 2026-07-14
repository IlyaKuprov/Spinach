% 2D microfluidic concentration plotting function. Uses mesh 
% tessellation information to plot concentrations as vertical
% bars. This function should be called after mesh_plot() has
% drawn the mesh. Syntax:
%
%              conc_plot(spin_system,conc,obs)
%
% Parameters:
%
%    spin_system - Spinach spin system object containing
%                  mesh and tessellation information
%
%    conc        - concentrations as a column vector with
%                  the same number of elements as the num-
%                  ber of Voronoi cells; these will deter-
%                  mine bar heights
%
%    obs         - up to three observables as columns of
%                  a matrix with the same number of rows
%                  as conc; these will be normalised and
%                  mapped into HSV colour space for each
%                  Voronoi cell bar. Options:
%
%                   one column:    [xy_phases]
%
%                   two columns:   [xy_phases xy_amps]
%
%                   three columns: [xy_phases xy_amps z]
%
% Outputs:
%
%    the function updates a figure created by mesh_plot()
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=conc_plot.m>

function conc_plot(spin_system,conc,obs)

% Check consistency
if exist('obs','var')
    grumble(spin_system,conc,obs);
else
    grumble(spin_system,conc);
end

% Decide the colours
if ~exist('obs','var')
    
    % Neutral grey if no observables supplied
    RGB=0.5*ones(spin_system.mesh.vor.ncells,3);

elseif size(obs,2)==1

    % One observable: assume phase and map into middle hues
    RGB=hsv2rgb(wrapTo2Pi(obs)/(2*pi),...
                0.75*ones(size(spin_system.mesh.vor.cells,1),1),...
                0.50*ones(size(spin_system.mesh.vor.cells,1),1));

elseif size(obs,2)==2

    % Two observables: assume phase + amp and map into HS
    RGB=hsv2rgb(wrapTo2Pi(obs(:,1))/(2*pi),...
                obs(:,2)/max(obs(:,2)),...
                0.50*ones(size(spin_system.mesh.vor.cells,1),1));

elseif size(obs,2)==3

    % Three observables: assume phase + amp + Z and map into HSV
    RGB=hsv2rgb(wrapTo2Pi(obs(:,1))/(2*pi),...
                obs(:,2)/max(obs(:,2)),...
                (obs(:,3)+min(obs(:,3)))/(max(obs(:,3))-min(obs(:,3))));

else

    % Complain and bomb out
    error('incorrect number of colour mapped observables.');

end

% Find cells with significant bar heights
active_cells=find(abs(conc)>1e-3*diff(spin_system.mesh.zext));
cell_sizes=cellfun(@numel,spin_system.mesh.vor.cells(active_cells));
total_vertices=sum(cell_sizes);

% Preallocate cap geometry and colours
V=zeros(total_vertices,3);
F=nan(numel(active_cells),spin_system.mesh.vor.max_cell_size+1);
FRGB=zeros(numel(active_cells),3);
vertex_offset=0;

% Build lids and bottoms
for m=1:numel(active_cells)

    % Get the vertices of the Voronoi cell
    n=active_cells(m); nvert=cell_sizes(m);
    vor_cell_x=spin_system.mesh.vor.vertices(spin_system.mesh.vor.cells{n},1);
    vor_cell_y=spin_system.mesh.vor.vertices(spin_system.mesh.vor.cells{n},2);
    vor_cell_z=conc(n)*ones(size(vor_cell_x));
    vertex_range=vertex_offset+(1:nvert);
    V(vertex_range,:)=[vor_cell_x(:) vor_cell_y(:) vor_cell_z(:)];

    % Build the face connectivity index
    face=[1:nvert 1]+vertex_offset;
    F(m,1:numel(face))=face;

    % Add the colour spec and advance the vertex offset
    FRGB(m,:)=RGB(n,:);
    vertex_offset=vertex_offset+nvert;

end

% Draw lids and bottoms of the bars as a multifaceted patch
patch('Faces',F,'Vertices',V,'FaceColor','flat',...
      'FaceVertexCData',FRGB,'EdgeColor','black',...
      'LineWidth',0.25,'LineJoin','round'); V(:,3)=0;
patch('Faces',F,'Vertices',V,'FaceColor','flat',...
      'FaceVertexCData',FRGB,'EdgeColor','black',...
      'LineWidth',0.125,'LineJoin','round');

% Preallocate side geometry and colours
V=zeros(2*total_vertices,3);
F=zeros(total_vertices,5);
FRGB=zeros(total_vertices,3);
vertex_offset=0; face_offset=0;

% Build sides
for m=1:numel(active_cells)

    % Get the vertices of the Voronoi cell
    n=active_cells(m); nvert=cell_sizes(m);
    vor_cell_x=spin_system.mesh.vor.vertices(spin_system.mesh.vor.cells{n},1);
    vor_cell_y=spin_system.mesh.vor.vertices(spin_system.mesh.vor.cells{n},2);
    vor_cell_z=conc(n)*ones(size(vor_cell_x));
    top_range=vertex_offset+(1:nvert);
    bottom_range=vertex_offset+nvert+(1:nvert);
    V(top_range,:)=[vor_cell_x(:) vor_cell_y(:) vor_cell_z(:)];
    V(bottom_range,:)=[vor_cell_x(:) vor_cell_y(:) zeros(nvert,1)];

    % Build all walls for this cell
    local_idx=(1:nvert)';
    next_idx=[(2:nvert)';1];
    face_range=face_offset+(1:nvert);
    F(face_range,:)=[local_idx next_idx next_idx+nvert ...
                     local_idx+nvert local_idx]+vertex_offset;
    FRGB(face_range,:)=repmat(RGB(n,:),nvert,1);

    % Advance the vertex and face offsets
    vertex_offset=vertex_offset+2*nvert;
    face_offset=face_offset+nvert;

end

% Draw the sides of the bars as a multifaceted patch
patch('Faces',F,'Vertices',V,'FaceColor','flat',...
      'FaceVertexCData',FRGB,'EdgeColor','black',...
      'LineWidth',0.25,'LineJoin','round');

end

% Consistency enforcement
function grumble(spin_system,conc,obs)
if ~isfield(spin_system,'mesh')
    error('mesh information is missing from the spin_system structure.');
end
if ~isfield(spin_system.mesh,'vor')
    error('Voronoi tessellation information is missing from spin_system.mesh structure.');
end
if (~isfield(spin_system.mesh,'zext'))||(~isnumeric(spin_system.mesh.zext))||...
   (~isreal(spin_system.mesh.zext))||(numel(spin_system.mesh.zext)~=2)||...
   any(~isfinite(spin_system.mesh.zext))
    error('spin_system.mesh.zext must be a finite real two-element vector.');
end
if (~isfield(spin_system.mesh.vor,'ncells'))||(~isfield(spin_system.mesh.vor,'cells'))||...
   (~isfield(spin_system.mesh.vor,'vertices'))||(~isfield(spin_system.mesh.vor,'max_cell_size'))
    error('Voronoi tessellation information is incomplete.');
end
if (~isnumeric(conc))||(~isreal(conc))||(~iscolumn(conc))||...
   any(~isfinite(conc))||(numel(conc)~=spin_system.mesh.vor.ncells)
    error('conc must be a finite real column vector with one value per Voronoi cell.');
end
if nargin>2
    if (~isnumeric(obs))||(~isreal(obs))||...
       any(~isfinite(obs(:)))||(size(obs,1)~=spin_system.mesh.vor.ncells)||...
       (size(obs,2)<1)||(size(obs,2)>3)
        error('obs must have one, two, or three real columns and one row per Voronoi cell.');
    end
end
end

% Their god is a corpse nailed to a tree.
%
% Vikings, about Christians,
% in The Northman (2022)
