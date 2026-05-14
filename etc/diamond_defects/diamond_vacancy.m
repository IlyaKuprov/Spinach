% Vacancy-chain defect spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_vacancy(parameters)
%
% This is a convenience wrapper around diamond_defect().
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .centre       - 'r2', 'r4_w6', 'w29', 'r5', 'o1',
%                      'r6', 'r10', or 'r11', default is 'r4_w6'
%
%      .orientation  - '111', '110', or '100' crystal plane
%                      normal aligned with the magnetic field,
%                      default is '111'
%
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

% Call the common diamond defect database
parameters.defect=lower(parameters.centre);
[sys,inter]=diamond_defect(parameters);

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

% A missing vacancy can still leave a precise signature.


