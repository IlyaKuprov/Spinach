% NV- center ground state spin system for diamond (15N). Syntax:
%
%                    [sys,inter]=diamond_nvm_gs(parameters)
%
% Magnetic parameters are according to 
%       S. Felton et al., Phys. Rev. B 79 (2009) 075203, 
% spin-lattice relaxation according to eq. 1 in 
%       A. Jarmola et al., PRL 108 (2012) 197601,
% electron spin T2 times assumed to be equal to T2=0.5*T1 as in CPMG
% experiments measured in 
%       N. Bar-Gill et al., Nat. Comm. 4 (2013) 1743
% in the strong spin-phonon coupling regime. 15N isotope is assumed, 
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

function [sys,inter]=diamond_nvm_gs_15n(parameters)

% Check consistency.
grumble(parameters);

% Set default concentration.
if~isfield(parameters,'concentration')
    parameters.concentration=0.001;
end

% Set default temperature.
if~isfield(parameters,'temperature')
    parameters.temperature=296;
end

% Set default orientation.
if~isfield(parameters,'orientation')
    parameters.orientation='111';
end

% Define spin system isotopes.
sys.isotopes={'E3','15N'};

% Compute electron relaxation rates.
% A. Jarmola et al., PRL 108 (2012) 197601, eq. (1)
a1=0.8*parameters.concentration;    % s^-1*ppm^-1
a2=2.1e3;                           % s^-1
a3=2.2e-11;                         % K^-5*s^-1
delta=847.1303254502;               % K
r1e=a1+a2/(exp(delta/parameters.temperature)-1)+a3*(parameters.temperature)^5;

% Set long T2, as in dynamical decoupling experiments
r2e=2.0*r1e;
% Alternatively, can be set as Tm~600us, as reported for natural abundance
% 13C diamonds in Stanwix et al., PRB 2010

% Compute nuclear relaxation rates.
% Let's take that John Morton's assumption for 31P in silicon is valid,
% T2n~2*T1e.
r2n=0.5*r1e;

% order of magnitude, consistent with Pfender et al., Nat.Comm. 8 (2017) 834
% and Soshenko et al. Quant. Electronics 51 (2021) 1144; although these 
% values relate to room temperature only.
r1n=1e-2;

% Set relaxation parameters.
inter.relaxation={'t1_t2'};
inter.r1_rates={r1e r1n};
inter.r2_rates={r2e r2n};
inter.equilibrium='zero';
inter.rlx_keep='diagonal';

% Define electron g-tensor.
g_e=diag([2.0031 2.0031 2.0029]);

% Define zero-field splitting tensor.
d_zfs=2872e6;
zfs_mat=diag([-d_zfs/3 -d_zfs/3 2*d_zfs/3]);

% Define hyperfine tensor.
a_hfc=diag([3.65 3.65 3.03])*1e6;

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
zfs_mat=rot_mat*zfs_mat*rot_mat';
a_hfc=rot_mat*a_hfc*rot_mat';

% Populate interaction specification.
inter.zeeman.matrix=cell(1,2);
inter.zeeman.matrix{1}=g_e;
inter.zeeman.matrix{2}=zeros(3);
inter.coupling.matrix=cell(2,2);
inter.coupling.matrix{1,1}=zfs_mat;
inter.coupling.matrix{1,2}=a_hfc;
inter.coupling.matrix{2,2}=zeros(3);

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

% Check temperature if present.
if isfield(parameters,'temperature')
    if(~isnumeric(parameters.temperature))||(~isreal(parameters.temperature))||...
       (~isscalar(parameters.temperature))||(parameters.temperature<=0)
        error('parameters.temperature must be a positive real scalar.');
    end
end

% Check concentration if present.
if isfield(parameters,'concentration')
    if(~isnumeric(parameters.concentration))||(~isreal(parameters.concentration))||...
       (~isscalar(parameters.concentration))||(parameters.concentration<=0)
        error('parameters.concentration must be a positive real scalar.');
    end
end

% Check orientation if present.
if isfield(parameters,'orientation')
    if(~ischar(parameters.orientation))
        error('parameters.orientation must be a character string.');
    end
end

end
