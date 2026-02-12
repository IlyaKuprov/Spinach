% P1 center spin system for diamond, 15N enriched. Syntax:
%
%                    [sys,inter]=diamond_p1_15n(parameters)
% Magnetic resonance parameters from Nir-Orad et al. PCCP 2024
%
%
% Parameters:
%
%   parameters.temperature  - temperature, K (default 296)
%
%   parameters.concentration - defect concentration, ppm (default 0.001)
%
%   parameters.orientation  - '111', '110', or '100' crystal plane
%                              normal aligned with the magnetic field
%
% Outputs:
%
%   sys   - Spinach system specification structure
%
%   inter - Spinach interaction specification structure
%
% ilya.kuprov@weizmann.ac.il
% alexey.bogdanov@weizmann.ac.il

function [sys,inter]=diamond_p1_15n(parameters)

% Check consistency.
grumble(parameters);

% Set default orientation.
if~isfield(parameters,'orientation')
    parameters.orientation='111';
end

% Define spin system isotopes.
sys.isotopes={'E','15N'};

% Compute electron relaxation rates.
% Carroll et al. JMR 2021
r1e=1/140e-6;   % T1~140 us
r2e=1/1.9e-6;   % Tm~1.9 us  - values for the room temp on the central line

% Compute nuclear relaxation rates.
r2n=0.5*r1e;    % T2n defined by T1e
r1n=1e-2;       % just a guess, ~100s

% Set relaxation parameters.
inter.relaxation={'t1_t2'};
inter.r1_rates={r1e r1n};
inter.r2_rates={r2e r2n};
inter.equilibrium='zero';
inter.rlx_keep='diagonal';

% Define electron g-tensor.
% Nir-Orad et al. PCCP 2024
g_e=diag([2.00220 2.00220 2.00218]);

% Define hyperfine tensor, converted from 14N, given 
% gamma(15N)/gamma(14N)=-1.4027, and no direct measurements reported yet.
% Nir-Orad et al. PCCP 2024
a_hfc=(-1.4027)*diag([81.3 81.3 114.0])*1e6;
% Smith et al. Phys. Rev. 115, 1546 (1959)
% a_hfc=(-1.4027)*diag([81.9 81.9 114.3])*1e6*;

% Build rotation matrix for orientation.
if strcmp(parameters.orientation,'100')
    rot_mat=rotmat_align([1 1 1],[1 0 0]);
elseif strcmp(parameters.orientation,'110')
    rot_mat=rotmat_align([1 1 1],[1 1 0]);
elseif strcmp(parameters.orientation,'111')
    rot_mat=eye(3);
else
    error('parameters.orientation must be ''111'', ''110'', or ''100''.');
end

% Rotate interaction tensors.
g_e=rot_mat*g_e*rot_mat';
a_hfc=rot_mat*a_hfc*rot_mat';

% Populate interaction specification.
inter.zeeman.matrix=cell(1,2);
inter.zeeman.matrix{1}=g_e;
inter.zeeman.matrix{2}=zeros(3);
inter.coupling.matrix=cell(2,2);
inter.coupling.matrix{1,1}=zeros(3);
inter.coupling.matrix{1,2}=a_hfc;

end

% Rotation matrix aligning two vectors.
function rot_mat=rotmat_align(v_from,v_to)

% Normalize the input vectors.
v_from=v_from/norm(v_from);
v_to=v_to/norm(v_to);

% Compute the rotation axis and angle.
axis_vec=cross(v_from,v_to);
sin_ang=norm(axis_vec);
cos_ang=dot(v_from,v_to);

% Handle parallel vectors.
if sin_ang<1e-12
    if cos_ang>0
        rot_mat=eye(3);
        return
    else
        axis_vec=null(v_from.');
        axis_vec=axis_vec(:,1);
        sin_ang=0;
        cos_ang=-1;
    end
end

% Normalize the rotation axis.
axis_vec=axis_vec/sin_ang;

% Build the cross-product matrix.
skew=[0 -axis_vec(3) axis_vec(2);...
      axis_vec(3) 0 -axis_vec(1);...
     -axis_vec(2) axis_vec(1) 0];

% Build the rotation matrix.
rot_mat=eye(3)+sin_ang*skew+(1-cos_ang)*(skew*skew);

end

% Consistency enforcement.
function grumble(parameters)

% Check the input type.
if(~isstruct(parameters))
    error('parameters must be a structure.');
end

% Check orientation if present.
if isfield(parameters,'orientation')
    if(~ischar(parameters.orientation))
        error('parameters.orientation must be a character string.');
    end
end

end
