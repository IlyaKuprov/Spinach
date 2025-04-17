% Hydrodynamic flow generator on a mesh. Builds diffusion and 
% flow generator using the mesh parameters in the spin_system
% object. Syntax:
%
%              F=flow_gen(spin_system,parameters)
%
% Parameters:
%
%    spin_system - Spinach system descriptor object 
%                  containing mesh subfields produ-
%                  ced by COMSOL import functions
%
%    parameters.diff - diffusion coefficient, m^2/s
%
% Outputs:
%
%    F - spatial motion generator matrix with the
%        dimension equal to the number of Voronoi
%        cells of the mesh
%
% a.acharya@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=flow_gen.m>

function F=flow_gen(spin_system,parameters)

% Default is zero diffusion coefficient
if ~isfield(parameters,'diff'), parameters.diff=0; end

% Check consistency
grumble(spin_system,parameters);

% Substructure pull for parfor
mesh=spin_system.mesh; ncells=mesh.vor.ncells;

% Take Voronoi cell areas into account
A_forw=spdiags(mesh.vor.weights,0,ncells,ncells);
A_back=spdiags(1./mesh.vor.weights,0,ncells,ncells);

% Build flow index
F=cell(ncells,1);
parfor k=1:ncells %#ok<*PFBNS>

    % Pull out cell area
    A_k=mesh.vor.weights(k); 
    
    % Find the parent vertex
    parent_vertex=mesh.idx.active(k);

    % Find triangles containing parent vertex
    nearby_triangles=find(sum(mesh.idx.triangles==parent_vertex,2));
    nearby_triangles=mesh.idx.triangles(nearby_triangles,:);

    % Collect all nearby vertices
    nearby_vertices=setdiff(nearby_triangles(:),parent_vertex);

    % Find the Voronoi cells they are parenting
    nearby_cells=find(ismember(mesh.idx.active,nearby_vertices))';
    
    % Start local flow table
    F_local=zeros(0,3);

    % Loop over nearby cells
    for m=nearby_cells

        % See if they share a boundary
        shared_pts=intersect(mesh.vor.cells{k},...
                             mesh.vor.cells{m});
        if numel(shared_pts)==2

            % Determine boundary length
            b_km=norm(mesh.vor.vertices(shared_pts(1),:)-...       % #NORMOK
                      mesh.vor.vertices(shared_pts(2),:),2);       % m

            % Distance vector between Voronoi cell centres
            r_km=[mesh.x(mesh.idx.active(m)); mesh.y(mesh.idx.active(m))]-...
                 [mesh.x(mesh.idx.active(k)); mesh.y(mesh.idx.active(k))]; 

            % Velocity in the current and the adjacent cell
            v_m=[mesh.u(mesh.idx.active(m)); mesh.v(mesh.idx.active(m))];
            v_k=[mesh.u(mesh.idx.active(k)); mesh.v(mesh.idx.active(k))];

            % Contribution to the off-diag part
            F_km=+(1/A_k)*(b_km/norm(r_km,2))*parameters.diff ...
                 -(1/A_k)*(b_km/norm(r_km,2))*dot((v_m+v_k)/2,r_km)/2;

            % Add the terms to the generator
            if F_km>0
                F_local=[F_local; k m  F_km];
            else
                F_local=[F_local; m k -F_km];
            end

        end

    end

    % Add to global
    F{k}=F_local;

end

% Build and balance the flow generator
F=cell2mat(F); F=sparse(F(:,1),F(:,2),F(:,3),ncells,ncells);
F=F-spdiags(sum(F,1)',0,ncells,ncells); F=A_back*F*A_forw;  

end

% Consistency enforcement
function grumble(spin_system,parameters)
if ~isfield(spin_system,'mesh')
    error('mesh information is missing from the spin_system structure.');
end
if ~isfield(spin_system.mesh,'idx')
    error('indexing information is missing from spin_system.mesh structure.');
end
if ~isfield(spin_system.mesh,'vor')
    error('Voronoi tessellation information is missing from spin_system.mesh structure.');
end
if ~isfield(parameters,'diff')
    error('diffusion coefficient must be specified in parameters.diff variable.');
end
end

% A disbelieving production team on a BBC Radio 4 arts programme
% once heard one of their contributors use a c-word during enthu-
% siastic praise for a video installation which -- in their view
% -- symbolised a part of the female anatomy. Nervously checking
% the call log later, the team were amazed to find a total lack 
% of response. Until, that is, they reached the listener who had 
% been "absolutely disgusted to hear someone on the programme 
% use a split infinitive".
%
% https://www.spectator.co.uk/article/the-art-of-swearing/

