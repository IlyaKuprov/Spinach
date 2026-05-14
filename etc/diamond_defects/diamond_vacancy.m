% Vacancy-family spin system dispatcher for diamond. Syntax:
%
%          [sys,inter]=diamond_vacancy(parameters)
%
% This compatibility function dispatches to specific vacancy and
% vacancy-chain constructors. New code should call those constructors
% directly.
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .centre       - 'r2', 'r4_w6', 'w6', 'r4', 'w29', 'r5',
%                      'o1', 'r6', 'r10', or 'r11', default is 'r4_w6'
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field, default is '111'
%      .d_sign       - sign of D for R2, default is +1
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_vacancy.m>

function [sys,inter]=diamond_vacancy(parameters)

% Set default input
if nargin==0
    parameters=struct();
end

% Check consistency
grumble(parameters);

% Set default centre
if ~isfield(parameters,'centre')
    parameters.centre='r4_w6';
end

% Dispatch to the requested vacancy constructor
centre=lower(parameters.centre);
switch centre
    case 'r2'
        [sys,inter]=diamond_r2(parameters);
    case {'r4_w6','w6','r4'}
        [sys,inter]=diamond_r4_w6(parameters);
    case 'w29'
        [sys,inter]=diamond_w29(parameters);
    case 'r5'
        [sys,inter]=diamond_r5(parameters);
    case 'o1'
        [sys,inter]=diamond_o1(parameters);
    case 'r6'
        [sys,inter]=diamond_r6(parameters);
    case 'r10'
        [sys,inter]=diamond_r10(parameters);
    case 'r11'
        [sys,inter]=diamond_r11(parameters);
    otherwise
        error('unknown vacancy centre.');
end

end

% Consistency enforcement
function grumble(parameters)
if(~isstruct(parameters))
    error('parameters must be a structure.');
end
if isfield(parameters,'centre')&&(~ischar(parameters.centre))
    error('parameters.centre must be a character string.');
end
end

% A vacancy dispatcher should not contain the vacancy catalogue.

