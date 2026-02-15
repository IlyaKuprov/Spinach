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
% The following fields exist the parameters.* structure:
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

% Electron relaxation rates
a1=0.8*parameters.concentration;    % s^-1*ppm^-1
a2=2.1e3;                           % s^-1
a3=2.2e-11;                         % K^-5*s^-1
delta=847.1303254502;               % K
r1e=a1+a2/(exp(delta/parameters.temperature)-1)+...
       a3*(parameters.temperature)^5;
r2e=2.0*r1e;

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
inter.zeeman.matrix{1}=R*diag([2.0031 2.0031 2.0029])*R';

% Nuclear shielding tensor is ignored
inter.zeeman.matrix{2}=zeros(3);

% Electron ZFS tensor
inter.coupling.matrix{1,1}=R*zfs2mat(2872e6,0,0,0,0)*R';

% HFC and NQI tensor
switch parameters.nitrogen

    case '14N'
        
        % Hyperfine and quadrupolar coupling
        inter.coupling.matrix{1,2}=R*diag([-2.70e6 -2.70e6 -2.14e6])*R';
        inter.coupling.matrix{2,2}=R*zfs2mat(-5.01e6,0,0,0,0)*R';

    case '15N'

        % Only hyperfine coupling, opposite sign
        inter.coupling.matrix{1,2}=R*diag([3.65e6 3.65e6 3.03e6])*R';
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

% Once, when sheltering under a tree during a storm near Lichfield,
% he was asked to marry a heavily pregnant bride to a rather guilty
% looking groom. Asked to provide evidence that he had performed
% the shotgun wedding, Swift found a piece of paper and wrote:
%
%   Under an oak, in stormy weather,
%   I joined this rogue and whore together;
%   And none but He who rules the thunder
%   Can put this rogue and whore asunder.
%
% Madeline Grant, about Jonathan Swift

