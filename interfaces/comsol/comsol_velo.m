% Imports ASCII 2D flow velocity files produced by COMSOL. Syntax:
%
%      spin_system=comsol_velo(spin_system,file_name)
%
% Parameters:
%
%    spin_system  - Spinach spin system object
%
%    file_name    - a character string
%
% Ouputs:
%
%    the following fields are added to spin_system object
%
%       mesh.u, mesh.v      - column vectors with velocities
%                             at each vertex of the mesh
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=comsol_velo.m>

function spin_system=comsol_velo(spin_system,file_name)

% Check consistency
grumble(spin_system,file_name);

% Open the file
fid=fopen(file_name,'r');

% Velocity readout count
while true
    next_line=fgetl(fid);
    if contains(next_line,'% Nodes:')
        nvel=textscan(next_line,'%s %s %f');
        nvel=nvel{3}; fgetl(fid); fgetl(fid);
        fgetl(fid); fgetl(fid); break
    end
end
report(spin_system,[num2str(nvel) ' velocity readouts']);

% Parse velocity readouts
X=nan(nvel,1); Y=nan(nvel,1);
U=nan(nvel,1); V=nan(nvel,1);
for n=1:nvel
    VS=textscan(fgetl(fid),'%f %f %f %f %f');
    X(n)=VS{1}; Y(n)=VS{2};
    U(n)=VS{4}; V(n)=VS{5}; 
end

% Close the file
fclose(fid);

% Check that vertex locations are the same
if (norm(spin_system.mesh.x-X,1)>1e-6)||...
   (norm(spin_system.mesh.y-Y,1)>1e-6)
    error('vertex locations are different in the velocity file.');
end

% Store velocities in the spin system object
spin_system.mesh.u=U; spin_system.mesh.v=V;

end

% Consistency enforcement
function grumble(spin_system,file_name)
if ~ischar(file_name)
    error('file_name must be a character string.');
end
if ~isfield(spin_system,'mesh')
    error('mesh information is missing from the spin_system structure.');
end
end

% Words have no power to impress the mind without
% the exquisite horror of their reality.
% 
% Edgar Allan Poe

