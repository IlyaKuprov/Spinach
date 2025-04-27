% Marks 2D microfluidic mesh vertices as inactive in hydrodyna-
% mic and diffusive transport processes. Syntax:
%
%                mesh=mesh_inact(mesh,vertex_list)
%
% Parameters:
%
%    mesh        - Spinach mesh object
%
%    vertex_list - row vector of integers specifying 
%                  the vertices to be inactivated
%
% Outputs:
%
%    mesh        - updated mesh object
%   
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=mesh_inact.m>

function mesh=mesh_inact(mesh,vertex_list)

% Check consistency
grumble(mesh,vertex_list);

% Update the active vertex list
mesh.idx.active=setdiff(mesh.idx.active,vertex_list);

% Zero out velocities and concentrations, if present
vertex_list=setdiff(1:numel(mesh.x),mesh.idx.active);
if isfield(mesh,'u'), mesh.u(vertex_list)=0; end
if isfield(mesh,'v'), mesh.v(vertex_list)=0; end
if isfield(mesh,'c'), mesh.c(vertex_list,:)=0; end

end

% Consistency enforcement
function grumble(mesh,vertex_list)
if ~isfield(mesh,'idx')
    error('indexing information is missing from mesh structure.');
end
if (~isnumeric(vertex_list))||(~isreal(vertex_list))||...
   (~isrow(vertex_list))||any(mod(vertex_list,1)~=0)||any(vertex_list<1)
   error('vertex_list must be a row vector of positive integers.');
end
if any(vertex_list>numel(mesh.x))
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

