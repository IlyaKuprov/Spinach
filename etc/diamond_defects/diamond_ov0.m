% Neutral oxygen-vacancy (OV0, WAR5) centre ground state spin system
% for diamond. Syntax:
%
%              [sys,inter]=diamond_ov0(parameters)
%
% Magnetic parameters from Tables 9-2 and 9-3 of:
%
%     B.L. Cann, Magnetic Resonance Studies of Point Defects in
%     Diamond, PhD thesis, University of Warwick (2009),
%     <https://wrap.warwick.ac.uk/id/eprint/3125>
%
% The zero-field splitting is confirmed at 4 K, and the defect is
% assigned to OV0, in:
%
%     S. Mukherjee et al., Phys. Rev. B 114, 074105 (2026),
%     <https://doi.org/10.1103/3dcd-mkcq>
%
% The centre has S=1 and C3v symmetry; the electron Zeeman and zero-
% field splitting tensors are axial about the trigonal axis, and the
% rhombicity is zero within the experimental error. Oxygen is left out
% of the spin system because 16O, which makes up the overwhelming ma-
% jority of the natural isotope mixture, has no nuclear spin.
%
% Four 13C hyperfine tensors are resolved in the experimental spectra,
% labelled a, g, l, and c after the corresponding shells of the NV(-)
% centre; the thesis assigns 3, 6, 3, and 3 symmetry-equivalent carbon
% atoms to them respectively. One carbon per shell is returned here;
% the tensors of the remaining atoms of each shell are obtained by the
% C3v operations around the trigonal axis. For the a and g shells the
% symmetry-relaxed fit is used, which the thesis reports as the better
% of the two.
%
% Parameters:
%
%   a structure (parameters.*) with the following field:
%
%      .orientation   - '111', '110', or '100' crystal plane normal
%                       aligned with the magnetic field
%
% Outputs:
%
%   sys   - Spinach system specification structure
%
%   inter - Spinach interaction specification structure
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=diamond_ov0.m>

function [sys,inter]=diamond_ov0(parameters)

% Check input count
if nargin~=1
    error('exactly one input argument is required.');
end

% Check consistency
grumble(parameters);

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
        error('unknown orientation specification.');

end

% Define spin system isotopes and labels
sys.isotopes={'E3','13C','13C','13C','13C'};
sys.labels={'OV0','C_a','C_g','C_l','C_c'};
inter.zeeman.matrix=cell(1,numel(sys.isotopes));
inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes));

% Electron g-tensor
inter.zeeman.matrix{1}=C*frame*diag([2.0025 2.0025 2.0029])*frame'*C';

% Electron ZFS tensor
inter.coupling.matrix{1,1}=C*frame*zfs2mat(2888e6,0,0,0,0)*frame'*C';

% 13C hyperfine eigenvalues, MHz
hfc_vals=[197.4 117.3 118.2;
          17.5  11.7  13.0;
          12.6  8.5   8.5;
          7.4   4.3   4.3];

% 13C hyperfine principal-axis polar angles from [001], degrees
hfc_theta=[53.86 143.86 90.00;
           60.40 150.40 90.00;
           54.70 144.74 90.00;
           54.70 144.74 90.00];

% 13C hyperfine principal-axis azimuths from [100], degrees
hfc_phi=[225.26 225.00 135.00;
         225.26 225.00 135.00;
         225.26 225.00 135.00;
         225.26 225.00 135.00];

% Build and rotate the 13C hyperfine tensors
for n=1:size(hfc_vals,1)

    % Convert principal-axis directions from polar angles
    axis_xyz=[sind(hfc_theta(n,:)).*cosd(hfc_phi(n,:));...
              sind(hfc_theta(n,:)).*sind(hfc_phi(n,:));...
              cosd(hfc_theta(n,:))];

    % Build the crystal-frame hyperfine tensor
    hfc_cub=axis_xyz*diag(1e6*hfc_vals(n,:))*axis_xyz';

    % Symmetrise and rotate into the requested orientation
    inter.coupling.matrix{1,n+1}=C*((hfc_cub+hfc_cub')/2)*C';

end

end

% Consistency enforcement
function grumble(parameters)
if(~isstruct(parameters))
    error('parameters must be a structure.');
end
if ~isfield(parameters,'orientation')
    error('parameters.orientation field is required.');
end
if(~ischar(parameters.orientation))
    error('parameters.orientation must be a character string.');
end
if ~ismember(parameters.orientation,{'111','110','100'})
    error('parameters.orientation must be ''111'', ''110'', or ''100''.');
end
end

% The best that most of us can hope to achieve in physics is simply
% to misunderstand at a deeper level.
%
% Wolfgang Pauli, attributed

