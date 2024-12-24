% ORCA cube file parser. Extracts the normalised probability density and
% the associated metric information from ORCA spin density in "3D simple
% format" (see ORCA manual). Syntax:
%
%            [density,ext,dx,dy,dz]=ocparse(filename,pad_factor)
%
% Parameters:
%
%    filename   - character string specifying the file to load
%
%    pad_factor - padding factor specifying how many multiples
%                 the array dimension in zeros to add on each 
%                 side of the cube
%
% Outputs:
%
%    density  - probability density cube with dimensions
%               ordered as [X Y Z]
%
%    ext      - grid extents in Angstrom, ordered as
%               [xmin xmax ymin ymax zmin zmax]
%
%    dx,dy,dz - grid steps in the three directions, Angstrom
%
% Ilya Kuprov (ilya.kuprov@weizmann.ac.uk)
% Elizaveta Suturina( e.suturina@soton.ac.uk)
% Petra Pikulova (484677@mail.muni.cz)
%
% <https://spindynamics.org/wiki/index.php?title=ocparse.m>

function [density,ext,dx,dy,dz]=ocparse(filename,pad_factor)

% Check consistency
grumble(filename,pad_factor);

% Inform the user
disp('Parsing ORCA cube...');

% Parse the file
A=importdata(filename,' ',4);

% Get the number of points along [x y z]
npts=str2num(A.textdata{2}); %#ok<ST2NM>

% Make the probability density cube
density=reshape(abs(A.data),npts(2),npts(3),npts(1));
density=permute(density,[3 2 1]);

% Get the corner coordinates 
corner_xyz=str2num(A.textdata{3}); %#ok<ST2NM>

% Get the grid spacing vector
dxdydz=str2num(A.textdata{4}); %#ok<ST2NM>
dx=dxdydz(1); dy=dxdydz(2); dz=dxdydz(3);

% Integrate and normalise probability density
total_prob=trapz(trapz(trapz(density)))*dx*dy*dz;
density=density/total_prob;

% Compute grid extents
ext=[corner_xyz(1) (corner_xyz(1)+(npts(1)-1)*dx)...
     corner_xyz(2) (corner_xyz(2)+(npts(2)-1)*dy)...
     corner_xyz(3) (corner_xyz(3)+(npts(3)-1)*dz)];
 
% Pad the density with zeros
density=padarray(density,pad_factor*size(density),0,'both');
ext=[ext(1)-abs(ext(2)-ext(1))*pad_factor ext(2)+abs(ext(2)-ext(1))*pad_factor...
     ext(3)-abs(ext(4)-ext(3))*pad_factor ext(4)+abs(ext(4)-ext(3))*pad_factor...
     ext(5)-abs(ext(6)-ext(5))*pad_factor ext(6)+abs(ext(6)-ext(5))*pad_factor];

end

% Consistency enforcement
function grumble(filename,pad_factor)
if ~exist(filename,'file')
    error('the file specified as not found.');
end
if (~isnumeric(pad_factor))||(~isreal(pad_factor))||...
   (~isscalar(pad_factor))||(pad_factor<0)||...
   (mod(pad_factor,1)~=0)
    error('pad_factor must be a non-negative real integer.');
end
end

% There is a cult of ignorance in the United States, and there
% has always been. The strain of anti-intellectualism has been
% a constant thread winding its way through our political and
% cultural life, nurtured by the false notion that democracy
% means that 'my ignorance is just as good as your knowledge'. 
%
% Isaac Asimov

