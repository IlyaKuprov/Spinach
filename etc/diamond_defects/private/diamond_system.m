% Spinach system builder for diamond defects. Syntax:
%
%          [sys,inter]=diamond_system(electron,gmat,zfs,nuclei,orientation)
%
% This internal helper builds Spinach sys and inter structures from
% electronic, nuclear, and tensor data supplied by diamond-defect
% constructors.
%
% Parameters:
%
%    electron     - Spinach electron isotope label
%
%    gmat         - electron g tensor in the diamond crystal frame
%
%    zfs          - electron zero-field tensor in Hz, or []
%
%    nuclei       - cell array of isotope, hyperfine, and quadrupole data
%
%    orientation  - '111', '110', or '100' crystal plane normal aligned
%                   with the magnetic field
%
% Outputs:
%
%    sys          - Spinach system specification structure
%
%    inter        - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_system.m>

function [sys,inter]=diamond_system(electron,gmat,zfs,nuclei,orientation)

% Set default orientation
if nargin<5
    orientation='111';
end

% Check consistency
grumble(electron,gmat,zfs,nuclei,orientation);

% Set crystal-to-laboratory rotation
C=diamond_orient(orientation);

% Initialise the electron spin
sys.isotopes={electron};
inter.zeeman.matrix{1}=C*gmat*C';

% Set electron zero-field splitting
if ~isempty(zfs)
    inter.coupling.matrix{1,1}=diamond_traceless(C*zfs*C');
else
    inter.coupling.matrix{1,1}=[];
end

% Add all nuclei
for n=1:numel(nuclei)
    sys.isotopes{n+1}=nuclei{n}.iso;
    inter.zeeman.matrix{n+1}=zeros(3);
    inter.coupling.matrix{1,n+1}=C*nuclei{n}.A*C';
    if ~isempty(nuclei{n}.Q)
        inter.coupling.matrix{n+1,n+1}=diamond_traceless(C*nuclei{n}.Q*C');
    else
        inter.coupling.matrix{n+1,n+1}=[];
    end
end

end

% Consistency enforcement
function grumble(electron,gmat,zfs,nuclei,orientation)
if ~ischar(electron)
    error('electron must be a character string.');
end
if(~isnumeric(gmat)||~isreal(gmat)||~isequal(size(gmat),[3 3]))
    error('gmat must be a real 3x3 matrix.');
end
if(~isempty(zfs)&&(~isnumeric(zfs)||~isreal(zfs)||~isequal(size(zfs),[3 3])))
    error('zfs must be empty or a real 3x3 matrix.');
end
if ~iscell(nuclei)
    error('nuclei must be a cell array.');
end
if ~ischar(orientation)
    error('orientation must be a character string.');
end
end

% A catalogue should not know why it is being searched.

