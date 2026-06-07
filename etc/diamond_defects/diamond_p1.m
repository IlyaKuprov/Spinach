% P1 centre spin system for diamond. Syntax:
%
%           [sys,inter]=diamond_p1(parameters)
%
% Magnetic parameters from: Nir-Orad et al. PCCP 2024, and
% Smith et al. Phys. Rev. 115, 1546 (1959), relaxation ra-
% tes from Carroll et al. JMR 2021.
%
% Parameters:
%
%   a structure (parameters.*) with the following fields:
%
%      .orientation   - '111', '110', or '100' crystal 
%                        plane normal aligned with the
%                        magnetic field, def. is '111'
%
%      .nitrogen      - '14N' or '15N', default is '14N'
%
% Outputs:
%
%   sys   - Spinach system specification structure
%
%   inter - Spinach interaction specification structure
%
% alexey.bogdanov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=diamond_p1.m>

function [sys,inter]=diamond_p1(parameters)

% Check consistency
grumble(parameters);

% Set default parameters
if ~isfield(parameters,'orientation')
    parameters.orientation='111';
end
if ~isfield(parameters,'nitrogen')
    parameters.nitrogen='14N';
end

% Define spin system isotopes
sys.isotopes={'E',parameters.nitrogen};
inter.zeeman.matrix=cell(1,numel(sys.isotopes));
inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes));

% Set the trigonal principal-axis frame
frame=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
        1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
        0          2/sqrt(6)  1/sqrt(3)];

% Rotation matrix for orientation
switch parameters.orientation
    
    case '100'

        C=rotmat_align([1 0 0],[0 0 1]);

    case '110'

        C=rotmat_align([1 1 0],[0 0 1]);

    case '111'

        C=rotmat_align([1 1 1],[0 0 1]);

    otherwise

        % Complain and bomb out
        error('unknown oritentation specification');

end

% Electron g-tensor
inter.zeeman.matrix{1}=C*frame*diag([2.00220 2.00220 2.00218])*frame'*C';

% HFC and NQI tensor
switch parameters.nitrogen

    case '14N'
        
        % Hyperfine and quadrupolar coupling
        inter.coupling.matrix{1,2}=C*frame*diag([81.3e6 81.3e6 114.0e6])*frame'*C';
        inter.coupling.matrix{2,2}=C*frame*zfs2mat(-3.97e6,0,0,0,0)*frame'*C';

    case '15N'

        % Only hyperfine coupling, opposite sign
        inter.coupling.matrix{1,2}=C*frame*diag([-114.0e6 -114.0e6 -159.9e6])*frame'*C';

    otherwise

        % Complain and bomb out
        error('wrong nitrogen isotope.');

end

end

% Consistency enforcement
function grumble(parameters)
if(~isstruct(parameters))
    error('parameters must be a structure.');
end
if isfield(parameters,'orientation')
    if(~ischar(parameters.orientation))
        error('parameters.orientation must be a character string.');
    end
end
end
