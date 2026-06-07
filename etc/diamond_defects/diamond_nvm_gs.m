% NV centre ground state spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_nvm_gs(parameters)
%
% Magnetic parameters from:
%
%     S. Felton et al., Phys. Rev. B 79 (2009) 075203
%
% Spin-lattice relaxation according to Eq. 1 in:
%
%       A. Jarmola et al., PRL 108 (2012) 197601
%
% Electron spin T2 times assumed to be equal to T2=0.5*T1 as
% in CPMG experiments reported in:
%
%     N. Bar-Gill et al., Nat. Comm. 4 (2013) 1743
%
% in the strong spin-phonon coupling regime. Nuclear relaxa-
% tion rartes are intelligent guesses based on:
%
%        Pfender et al., Nat.Comm. 8 (2017) 834
%
%    Soshenko et al. Quant. Electronics 51 (2021) 1144
%
% Parameters: 
% 
%  the following is needed in the parameters.* structure:
%
%   .temperature   - temperature, K (default: 298)
%
%   .concentration - defect concentration in ppm for the 
%                    relaxation model, default is 0.001
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
% <https://spindynamics.org/wiki/index.php?title=diamond_nvm_gs.m>

function [sys,inter]=diamond_nvm_gs(parameters)

% Check consistency
grumble(parameters);

% Set default parameters
if ~isfield(parameters,'concentration')
    parameters.concentration=0.001;
end
if ~isfield(parameters,'temperature')
    parameters.temperature=298;
end
if ~isfield(parameters,'orientation')
    parameters.orientation='111';
end
if ~isfield(parameters,'nitrogen')
    parameters.nitrogen='14N';
end

% Define spin system isotopes
sys.isotopes={'E3',parameters.nitrogen};
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
inter.zeeman.matrix{1}=C*frame*diag([2.0031 2.0031 2.0029])*frame'*C';

% Electron ZFS tensor
inter.coupling.matrix{1,1}=C*frame*zfs2mat(2872e6,0,0,0,0)*frame'*C';

% HFC and NQI tensor
switch parameters.nitrogen

    case '14N'
        
        % Hyperfine and quadrupolar coupling
        inter.coupling.matrix{1,2}=C*frame*diag([-2.70e6 -2.70e6 -2.14e6])*frame'*C';
        inter.coupling.matrix{2,2}=C*frame*zfs2mat(-5.01e6,0,0,0,0)*frame'*C';

    case '15N'

        % Only hyperfine coupling, opposite sign
        inter.coupling.matrix{1,2}=C*frame*diag([3.65e6 3.65e6 3.03e6])*frame'*C';

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
if isfield(parameters,'temperature')
    if(~isnumeric(parameters.temperature))||(~isreal(parameters.temperature))||...
       (~isscalar(parameters.temperature))||(parameters.temperature<=0)
        error('parameters.temperature must be a positive real scalar.');
    end
end
if isfield(parameters,'concentration')
    if(~isnumeric(parameters.concentration))||(~isreal(parameters.concentration))||...
       (~isscalar(parameters.concentration))||(parameters.concentration<=0)
        error('parameters.concentration must be a positive real scalar.');
    end
end
if isfield(parameters,'orientation')
    if(~ischar(parameters.orientation))
        error('parameters.orientation must be a character string.');
    end
end
end

