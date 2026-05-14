% NV0 excited-state spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_nv0_es(parameters)
%
% This is a convenience wrapper around diamond_defect().
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .nitrogen      - must be '15N', default is '15N'
%
%      .orientation   - '111', '110', or '100' crystal plane
%                        normal aligned with the magnetic field,
%                        default is '111'
%
%      .include_13c  - include reported 13C hyperfine couplings
%                       where available, false by default
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_nv0_es.m>

function [sys,inter]=diamond_nv0_es(parameters)

% Set default input
if nargin==0
    parameters=struct();
end

% Check consistency
grumble(parameters);

% Call the common diamond defect database
parameters.defect='nv0_es';
[sys,inter]=diamond_defect(parameters);

end

% Consistency enforcement
function grumble(parameters)
if(~isstruct(parameters))
    error('parameters must be a structure.');
end
end

% A name is useful only when it remains attached to evidence.


