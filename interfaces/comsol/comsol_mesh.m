% Imports ASCII 2D mesh files produced by COMSOL. Syntax:
%
%      spin_system=comsol_mesh(spin_system,file_name)
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
%       mesh.x, mesh.y      - column vectors with vertex
%                             coordinates
%
%       mesh.idx.edges      - two-column array of integers
%                             containing edge index
%
%       mesh.idx.triangles  - three-column array of integers
%                             containing triangle index
%
%       mesh.idx.rectangles - four-column array of integers
%                             containing rectangle index
%
% ilya.kuprov@weizmann.ac.il
% a.acharya@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=comsol_mesh.m>

function spin_system=comsol_mesh(spin_system,file_name)

% Check consistency
grumble(file_name);

% Open the file
fid=fopen(file_name,'r');

% Point count
while true
    next_line=fgetl(fid);
    if contains(next_line,'number of mesh points')
        npts=textscan(next_line,'%f');
        npts=npts{1}; break
    end
end
report(spin_system,[num2str(npts) ' mesh vertices']);

% Scroll the file to point coordinates
while ~contains(fgetl(fid),'Mesh point coordinates'), end

% Read vertex coordinates
mesh.x=nan(npts,1); 
mesh.y=nan(npts,1);
for n=1:npts
    XY=textscan(fgetl(fid),'%f %f');
    mesh.x(n)=XY{1}; mesh.y(n)=XY{2};
end

% Scroll the file to edge specification
while ~contains(fgetl(fid),'edg # type name'), end

% Edge count
while true
    next_line=fgetl(fid);
    if contains(next_line,'number of elements')
        nedg=textscan(next_line,'%f');
        nedg=nedg{1}; fgetl(fid); break
    end
end
report(spin_system,[num2str(nedg) ' mesh edges']);

% Read edges
mesh.idx.edges=nan(nedg,2);
for n=1:nedg
    ES=textscan(fgetl(fid),'%f %f');
    mesh.idx.edges(n,1)=ES{1}+1; 
    mesh.idx.edges(n,2)=ES{2}+1;
end

% Scroll the file to triangle specification
while ~contains(fgetl(fid),'tri # type name'), end

% Triangle count
while true
    next_line=fgetl(fid);
    if contains(next_line,'number of elements')
        ntri=textscan(next_line,'%f');
        ntri=ntri{1}; fgetl(fid); break
    end
end
report(spin_system,[num2str(ntri) ' mesh triangles']);

% Read triangles
mesh.idx.triangles=nan(ntri,3);
for n=1:ntri
    TS=textscan(fgetl(fid),'%f %f %f');
    mesh.idx.triangles(n,1)=TS{1}+1;
    mesh.idx.triangles(n,2)=TS{2}+1; 
    mesh.idx.triangles(n,3)=TS{3}+1; 
end

% Scroll the file to rectangle specification
while ~contains(fgetl(fid),'quad # type name'), end

% Rectangle count
while true
    next_line=fgetl(fid);
    if contains(next_line,'number of elements')
        nrec=textscan(next_line,'%f');
        nrec=nrec{1}; fgetl(fid); break
    end
end
report(spin_system,[num2str(nrec) ' mesh rectangles']);

% Read rectangles 
mesh.idx.rectangles=nan(nrec,4);
for n=1:nrec
    RS=textscan(fgetl(fid),'%f %f %f %f');
    mesh.idx.rectangles(n,1)=RS{1}+1; 
    mesh.idx.rectangles(n,2)=RS{2}+1; 
    mesh.idx.rectangles(n,3)=RS{3}+1; 
    mesh.idx.rectangles(n,4)=RS{4}+1; 
end

% Add to the object
spin_system.mesh=mesh;

% Close the file
fclose(fid);

end

% Consistency enforcement
function grumble(file_name)
if ~ischar(file_name)
    error('file_name must be a character string.');
end
end

% Although Ray had a very much hands-off approach to men-
% toring, he occasionally would wander into the spectro-
% meter room and try to fine tune the shims on his newly
% acquired Varian XL-200 spectrometer, sometimes with di-
% sastrous consequences. Gareth Morris, who had returned
% as a postdoc to his group in 1979, solved this problem
% for us by installing a very large dial on the XL-200 
% front panel, not connected to anything on the inside,
% and marked in large lettering "Supervisor Fine Control".
%
% Adriaan Bax, about Ray Freeman

