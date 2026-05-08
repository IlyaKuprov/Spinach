%KEHL_ENDOR_AXIS ENDOR sweep axis for Kehl pulse sequences.
%
%   PARAMETERS=KEHL_ENDOR_AXIS(SPIN_SYSTEM,PARAMETERS,MODE) appends the
%   ENDOR Larmor frequencies and sweep-axis data used by Kehl sequence
%   kernels. MODE may be 'endor' for frequency-swept ENDOR, or 'time' for
%   the time-domain sequence.
%
%   Inputs:
%
%      SPIN_SYSTEM - Spinach spin system structure.
%      PARAMETERS  - Kehl ENDOR context parameters.
%      MODE        - 'endor' or 'time'.
%
%   Output:
%
%      PARAMETERS - input parameters with PARAMSENDOR appended.
%
%   February 2024 A. Kehl (akehl@gwdg.de)
%   Spinach architecture migration May 2026 Talos

function parameters=kehl_endor_axis(spin_system,parameters,mode)
if nargin<3
    mode='endor';
end

% Check consistency
grumble(spin_system,parameters,mode);

% Calculate nuclear Larmor frequencies
nuclei=parameters.endor_isotopes;
v_L=zeros(size(nuclei,2),1);
for n=1:size(nuclei,2)

    % Convert Spinach magnetogyric ratio to Hz/T
    [gamma,~]=spin(nuclei{n});
    v_L(n)=gamma*parameters.field_t/(2*pi);
end

% Build the requested ENDOR axis
if strcmp(mode,'time')

    % Time-domain x-axis uses microsecond-style legacy scaling
    range_EN=parameters.endor_range_hz*1e-12;
    res_EN=parameters.endor_res_hz*1e-12;
    Npts_EN=round(range_EN/res_EN);
    start_EN=parameters.time_start_s;
    step_EN=range_EN/(Npts_EN-1);
else

    % Frequency-domain x-axis uses the ENDOR RF sweep
    Npts_EN=round(parameters.endor_range_hz/parameters.endor_res_hz);
    if isfield(parameters,'rf_start_hz')
        start_EN=parameters.rf_start_hz;
    else
        start_EN=v_L(1)-parameters.endor_range_hz/2;
    end
    step_EN=parameters.endor_range_hz/(Npts_EN-1);
end

% Populate the Kehl ENDOR parameter map used by sequence kernels
x_coords=zeros(Npts_EN);
for n=1:Npts_EN
    x_coords(n)=start_EN+(n-1)*step_EN;
end
paramsENDOR=containers.Map;
paramsENDOR('v_L')=v_L;
paramsENDOR('start_EN')=start_EN;
paramsENDOR('step_EN')=step_EN;
paramsENDOR('range_EN')=parameters.endor_range_hz;
paramsENDOR('Npts_EN')=Npts_EN;
paramsENDOR('x_coords')=x_coords;
parameters.paramsENDOR=paramsENDOR;

end

% Consistency enforcement
function grumble(spin_system,parameters,mode)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if (~ischar(mode))&&(~isstring(mode))
    error('mode must be a character string or a string scalar.');
end
if ~ismember(char(mode),{'endor','time'})
    error('mode must be endor or time.');
end
end
