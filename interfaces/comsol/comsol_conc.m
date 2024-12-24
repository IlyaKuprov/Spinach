% Imports ASCII 2D concentration files produced by COMSOL. Syntax:
%
%      spin_system=comsol_conc(spin_system,file_name)
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
%       mesh.c    - stack of column vectors with 
%                   concentrations (in rows) at
%                   each vertex of the mesh
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=comsol_conc.m>

function spin_system=comsol_conc(spin_system,file_name)

% Check consistency
grumble(spin_system,file_name);

% Open the file
fid=fopen(file_name,'r');

% Concentration readout count
while true
    next_line=fgetl(fid);
    if contains(next_line,'% Nodes:')
        ncon=textscan(next_line,'%s %s %f');
        ncon=ncon{3}; fgetl(fid); fgetl(fid);
        fgetl(fid); fgetl(fid); break
    end
end
report(spin_system,[num2str(ncon) ' concentration readouts']);

% Read concentrations
X=nan(ncon,1); Y=nan(ncon,1); C=cell(ncon,1);
for n=1:ncon
    VS=textscan(fgetl(fid),'%f'); VS=VS{1};
    X(n)=VS(1); Y(n)=VS(2); C{n}=VS(4:end)'; 
end
spin_system.mesh.c=cell2mat(C);
report(spin_system,['for ' num2str(size(spin_system.mesh.c,2)) ' substances']);

% Close the file
fclose(fid);

% Check that vertex locations are the same
if (norm(spin_system.mesh.x-X,1)>1e-6)||...
   (norm(spin_system.mesh.y-Y,1)>1e-6)
    error('vertex locations are different in the concentrations file.');
end

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

% Feed a captive wolf as much as you like - he would still
% leg it into the forest at the first opportunity.
%
% A Russian saying

