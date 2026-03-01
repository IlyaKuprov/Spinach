% P1 centre spin system for diamond. Syntax:
%
%           [sys,inter]=diamond_p1(parameters)
%
% Magnetic parameters from: Nir-Orad et al. PCCP 2024, and
% Smith et al. Phys. Rev. 115, 1546 (1959), relaxation ra-
% tes from Carroll et al. JMR 2021.
%
% The following fields exist the parameters.* structure:
%
%   .orientation   - '111', '110', or '100' crystal 
%                     plane normal aligned with the
%                     magnetic field, def. is '111'
%
%   .nitrogen      - '14N' or '15N', default is '14N'
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

% Electron relaxation rates
% (room temp, central line)
r1e=1/140e-6;   % T1 ~ 140 µs
r2e=1/1.9e-6;   % T2 ~ 1.9 µs

% Nuclear relaxation rates
r1n=1e-2; r2n=0.5*r1e; 

% Relaxation parameters
inter.relaxation={'t1_t2'};
inter.r1_rates={r1e r1n};
inter.r2_rates={r2e r2n};
inter.equilibrium='zero';
inter.rlx_keep='diagonal';

% Rotation matrix for orientation
switch parameters.orientation
    
    case '100'

        R=rotmat_align([1 1 1],[1 0 0]);

    case '110'

        R=rotmat_align([1 1 1],[1 1 0]);

    case '111'

        R=eye(3);

    otherwise

        % Complain and bomb out
        error('unknown oritentation specification');

end

% Electron g-tensor
inter.zeeman.matrix{1}=R*diag([2.00220 2.00220 2.00218])*R';

% Nuclear shielding tensor is ignored
inter.zeeman.matrix{2}=zeros(3);

% HFC and NQI tensor
switch parameters.nitrogen

    case '14N'
        
        % Hyperfine and quadrupolar coupling
        inter.coupling.matrix{1,2}=R*diag([81.3e6 81.3e6 114.0e6])*R';
        inter.coupling.matrix{2,2}=R*zfs2mat(-3.97e6,0,0,0,0)*R';

    case '15N'

        % Only hyperfine coupling, opposite sign
        inter.coupling.matrix{1,2}=R*diag([-114.0e6 -114.0e6 -159.9e6])*R';
        inter.coupling.matrix{2,2}=[];

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

% One observes the survivors and learns from them.
%
% Frank Herbert, in the Dune series

