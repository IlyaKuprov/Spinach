% 2D microfluidic mesh plotting function. Draws the mesh, its Vo-
% ronoi tessellation, and a quiver plot of velocities. Syntax:
%
%           mesh_plot(spin_system,qscale,nodelabels)
%
% Parameters:
%
%    spin_system - Spinach spin system object containing
%                  mesh information
%
%    qscale      - scaling multiplier for the quiver plot
%                  of flow velocities, zero turns veloci-
%                  ty plotting off
%
%    nodelabels  - 1 causes vertex numbers to be displa-
%                  yed, 0 turns that off
%
% Outputs:
%
%    the function creates a figure
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=mesh_plot.m>

function mesh_plot(spin_system,qscale,nodelabels)

% Check consistency
grumble(spin_system,qscale,nodelabels);

% % Draw the edges
% patch(spin_system.mesh.plot.edg_a,...
%       spin_system.mesh.plot.edg_b,...
%       0,'EdgeColor',[0.5 0.5 0.5],'FaceColor','none',...
%         'LineWidth',0.125,'LineJoin','round');

% Draw the triangles
patch(spin_system.mesh.plot.tri_a,...
      spin_system.mesh.plot.tri_b,...
      0,'EdgeColor',[0.8 0.5 0.5],'FaceColor','none',...
        'LineWidth',0.125,'LineJoin','round'); 

% Set up the figure
hold on; box on; grid on; axis equal;

% Draw the rectangles
patch(spin_system.mesh.plot.rec_a,...
      spin_system.mesh.plot.rec_b,...
      0,'EdgeColor',[0.5 0.5 0.5],'FaceColor','none',...
        'LineWidth',0.125,'LineJoin','round');

% Draw Voronoi cells
patch(spin_system.mesh.plot.vor_a,...
      spin_system.mesh.plot.vor_b,...
      0,'EdgeColor',[0.5 0.8 0.5],'FaceColor','none',...
        'LineWidth',0.125,'LineJoin','round');

% Axis labels
kxlabel('X position, mm'); kylabel('Y position, mm');

% Draw quiver plot of velocities
if abs(qscale)>0
    quiver(spin_system.mesh.x,spin_system.mesh.y,...
           spin_system.mesh.u,spin_system.mesh.v,...
           qscale,'Color',[0 0.4470 0.7410]);
end

% Assign labels to mesh and wall points 
if nodelabels
    np_mesh=size(spin_system.mesh.x,1);
    plabels=arrayfun(@(n){sprintf('%d',n)},(1:np_mesh)');
    hpl=text(spin_system.mesh.x,...
             spin_system.mesh.y,plabels,'FontSize', 8);
    set(hpl,'Clipping','on');  
end

end

% Consistency enforcement
function grumble(spin_system,qscale,nodelabels)
if ~isfield(spin_system,'mesh')
    error('mesh information is missing from the spin_system structure.');
end
if ~isfield(spin_system.mesh,'plot')
    error('preplot information is missing from the spin_system structure.');
end
if (~isnumeric(qscale))||(~isreal(qscale))||...
   (~isscalar(qscale))||(qscale<0)
    error('qscale must be a non-negative real scalar.');
end
if (~isnumeric(nodelabels))||(~isreal(nodelabels))||...
   (~isscalar(nodelabels))||(~ismember(nodelabels,[0 1]))
    error('nodelabels must be 1 or 0.');
end
end

% Never attribute to malice that which is adequately 
% explained by stupidity.
%
% Hanlon's razor

