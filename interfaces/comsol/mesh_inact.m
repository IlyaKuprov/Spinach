% Marks 2D microfluidic mesh vertices as inactive in hydrodyna-
% mic and diffusive transport processes. Syntax:
%
%       spin_system=mesh_inact(spin_system,vertex_list)
%
% Parameters:
%
%    spin_system - Spinach data structure with a .mesh
%                  subfield present
%
%    vertex_list - row vector of integers specifying 
%                  the vertices to be inactivated
%
% Outputs:
%
%    spin_system - updated data structure
%   
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=mesh_inact.m>

function spin_system=mesh_inact(spin_system,vertex_list)

% Check consistency
grumble(spin_system,vertex_list);

% Update the active vertex list
spin_system.mesh.idx.active=setdiff(spin_system.mesh.idx.active,vertex_list);

% Zero out velocities and concentrations, if present
vertex_list=setdiff(1:numel(spin_system.mesh.x),spin_system.mesh.idx.active);
if isfield(spin_system.mesh,'u'), spin_system.mesh.u(vertex_list)=0; end
if isfield(spin_system.mesh,'v'), spin_system.mesh.v(vertex_list)=0; end
if isfield(spin_system.mesh,'c'), spin_system.mesh.c(vertex_list,:)=0; end

end

% Consistency enforcement
function grumble(spin_system,vertex_list)
if ~isfield(spin_system,'mesh')
    error('mesh information is missing from the spin_system structure.');
end
if ~isfield(spin_system.mesh,'idx')
    error('indexing information is missing from spin_system.mesh structure.');
end
if (~isnumeric(vertex_list))||(~isreal(vertex_list))||...
   (~isrow(vertex_list))||any(mod(vertex_list,1)~=0)||any(vertex_list<1)
   error('vertex_list must be a row vector of positive integers.');
end
if any(vertex_list>numel(spin_system.mesh.x))
    error('vertex index exceeds the number of vertices in the mesh.');
end
end

% The basic principle of the new education is to be that dunces and
% idlers must not be made to feel inferior to intelligent and indus-
% trious pupils. That would be "undemocratic". These differences be-
% tween pupils - for there are obviously and nakedly individual dif-
% ferences - must be disguised. This can be done at various levels.
% At universities, examinations must be framed so that nearly all the
% students get good marks. [...] But all the time there must be no
% faintest hint that they are inferior to the children who are at 
% work. Whatever nonsense they are engaged in must have - I believe
% the English already use the phrase - "parity of esteem".
%
% C.S. Lewis

