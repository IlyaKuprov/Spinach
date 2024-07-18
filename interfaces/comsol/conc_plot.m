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
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=conc_plot.m>

function conc_plot(spin_system,conc,obs)

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

% Check consistency
grumble(spin_system,conc);

% Loop over cells
V=zeros(0,3); patch_idx=0; FRGB=ones(0,3);
F=nan(0,spin_system.mesh.vor.max_cell_size+1);
for n=1:spin_system.mesh.vor.ncells

    % For significant bar heights
    if abs(conc(n))>1e-3*diff(spin_system.mesh.zext)

        % Get the vertices of the Voronoi cell
        vor_cell_x=spin_system.mesh.vor.vertices(spin_system.mesh.vor.cells{n},1);
        vor_cell_y=spin_system.mesh.vor.vertices(spin_system.mesh.vor.cells{n},2);
        vor_cell_z=conc(n)*ones(size(vor_cell_x));
        V=[V; vor_cell_x(:) vor_cell_y(:) vor_cell_z(:)]; %#ok<AGROW> 
        
        % Build the face connectivity index
        F_current=nan([1 spin_system.mesh.vor.max_cell_size+1]);
        new_face=[1:numel(vor_cell_z) 1]+patch_idx;
        F_current(1:numel(new_face))=new_face;
        F=[F; F_current]; patch_idx=new_face(end-1);      %#ok<AGROW>

        % Add the colour spec
        FRGB=[FRGB; RGB(n,:)];                            %#ok<AGROW>
        
    end
    
end

% Draw lids and bottoms of the bars as a multifaceted patch
patch('Faces',F,'Vertices',V,'FaceColor','flat',...
      'FaceVertexCData',FRGB,'EdgeColor','black',...
      'LineWidth',0.25,'LineJoin','round'); V(:,3)=0;
patch('Faces',F,'Vertices',V,'FaceColor','flat',...
      'FaceVertexCData',FRGB,'EdgeColor','black',...
      'LineWidth',0.125,'LineJoin','round');

% Loop over cells
V=zeros(0,3); patch_idx=0; 
F=nan(0,5); FRGB=ones(0,3);
for n=1:spin_system.mesh.vor.ncells

    % For significant bar heights
    if abs(conc(n))>1e-3*diff(spin_system.mesh.zext)

        % Get the vertices of the Voronoi cell
        vor_cell_x=spin_system.mesh.vor.vertices(spin_system.mesh.vor.cells{n},1);
        vor_cell_y=spin_system.mesh.vor.vertices(spin_system.mesh.vor.cells{n},2);
        vor_cell_z=conc(n)*ones(size(vor_cell_x));
        V=[V; vor_cell_x(:) vor_cell_y(:)   vor_cell_z(:); 
              vor_cell_x(:) vor_cell_y(:) 0*vor_cell_z(:)]; %#ok<AGROW> 
        nvert=numel(vor_cell_z);
        
        % Build the face connectivity index
        F_current=nan(nvert,5); FRGB_current=nan(nvert,3);
        for k=1:(nvert-1)

            % Build each wall
            F_current(k,:)=[k k+1 k+nvert+1 k+nvert k]+patch_idx;

            % Add the colour spec
            FRGB_current(k,:)=RGB(n,:);

        end
        F_current(nvert,:)=[nvert 1 nvert+1 2*nvert nvert]+patch_idx;
        FRGB_current(nvert,:)=RGB(n,:);
        
        % Add walls to the total
        patch_idx=patch_idx+2*nvert;
        F=[F; F_current]; FRGB=[FRGB; FRGB_current]; %#ok<AGROW>
        
    end
    
end

% Draw the sides of the bars as a multifaceted patch
patch('Faces',F,'Vertices',V,'FaceColor','flat',...
      'FaceVertexCData',FRGB,'EdgeColor','black',...
      'LineWidth',0.25,'LineJoin','round');

end

% Consistency enforcement
function grumble(spin_system,conc)
if ~isfield(spin_system,'mesh')
    error('mesh information is missing from the spin_system structure.');
end
if ~isfield(spin_system.mesh,'vor')
    error('Voronoi tessellation information is missing from spin_system.mesh structure.');
end
if (~isnumeric(conc))||(~isreal(conc))||(~iscolumn(conc))||...
   (numel(conc)~=spin_system.mesh.vor.ncells)
    error('conc must be a column vector with one real value per Voronoi cell.');
end
end

% Their god is a corpse nailed to a tree.
%
% Vikings, about Christians,
% in The Northman (2022)

