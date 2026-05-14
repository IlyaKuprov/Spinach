% Compatibility dispatcher for EPR-active defects in diamond. Syntax:
%
%          [sys,inter]=diamond_defect(parameters)
%
% This function preserves the original defect-label interface and
% dispatches to the specific diamond-defect constructors in this
% directory. New code should call the specific constructor directly.
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .defect       - defect label, for example 'nv0_es', 'n2vm',
%                      'war9', 'r2', 'siv0', 'gev0', 'w8', 'ne1',
%                      'n3', 'o4', 'ma1', or 'np9'
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field, default is '111'
%      .include_13c - include reported 13C hyperfine couplings where
%                     available, false by default
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_defect.m>

function [sys,inter]=diamond_defect(parameters)

% Check consistency
grumble(parameters);

% Select the specific constructor
defect=lower(parameters.defect);
switch defect
    case {'nv0_es','nv0'}
        [sys,inter]=diamond_nv0_es(parameters);
    case {'n2vm','n2v-'}
        [sys,inter]=diamond_n2vm(parameters);
    case {'war9','war10'}
        parameters.centre=defect;
        [sys,inter]=diamond_n_inter(parameters);
    case 'r2'
        [sys,inter]=diamond_r2(parameters);
    case {'r4_w6','w6','r4','w29','r5','o1','r6','r10','r11'}
        parameters.centre=defect;
        [sys,inter]=diamond_vacancy(parameters);
    case 'siv0'
        [sys,inter]=diamond_siv0(parameters);
    case 'gev0'
        [sys,inter]=diamond_gev0(parameters);
    case {'w8','ne1','ne2','ne3','ne4','ne5','ne8','ab1','ab2','ab3','ab4','ab5','nol1','nirim5'}
        parameters.centre=defect;
        [sys,inter]=diamond_ni(parameters);
    case {'n3','ok1'}
        parameters.centre=defect;
        [sys,inter]=diamond_ti(parameters);
    case {'o4','nlo2'}
        parameters.centre=defect;
        [sys,inter]=diamond_co(parameters);
    case {'ma1','np1','np2','np3','np4','np5','np6','np8','np9'}
        parameters.centre=defect;
        [sys,inter]=diamond_p(parameters);
    otherwise
        error('unknown diamond defect specification.');
end

end

% Consistency enforcement
function grumble(parameters)
if(~isstruct(parameters))
    error('parameters must be a structure.');
end
if ~isfield(parameters,'defect')
    error('parameters.defect must be specified.');
end
if(~ischar(parameters.defect))
    error('parameters.defect must be a character string.');
end
end

% A dispatcher should direct, not contain, the catalogue.

