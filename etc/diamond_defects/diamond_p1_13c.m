% P1 centre spin system with 13C neighbours in diamond. Syntax:
%
%           [sys,inter]=diamond_p1_13c(parameters)
%
% The electron and nitrogen parameters follow diamond_p1.m. The 13C
% hyperfine tensors and site assignments are from:
%
%   R C Barklie and J Guven, J. Phys. C: Solid State Phys. 14,
%   3621-3631 (1981), doi:10.1088/0022-3719/14/25/009
%
%   A Cox, M E Newton, and J M Baker, J. Phys.: Condens. Matter
%   6, 551-563 (1994), doi:10.1088/0953-8984/6/2/012
%
%   C V Peaker, M K Atumi, J P Goss, P R Briddon, A B Horsfall,
%   M J Rayson, and R Jones, Diamond Relat. Mater. 70, 118-123
%   (2016), doi:10.1016/j.diamond.2016.10.013
%
% Coordinates are representative site coordinates reconstructed from
% Peaker et al. Figure 1 and Tables 2-4. Their paper reports distances
% from the mid-point of the broken N-C bond in units of a0, but does
% not report relaxed Cartesian coordinates. The directions below use
% ideal diamond-lattice site directions, the radial distances use the
% Peaker et al. tabulated d/a0 values, and a0 is set to 3.567 Angstrom.
% The nitrogen atom is placed at the origin. The electron coordinate is
% left empty to avoid adding point-dipolar electron-nuclear couplings
% on top of the measured or calculated hyperfine tensors.
%
% The routine initialises the spin system with one nitrogen nucleus and
% eighteen 13C nuclei. Nuclei are labelled according to Peaker, 2016, and
% labels are listed in sys.labels
%
% Parameters:
%
%   a structure (parameters.*) with the following fields:
%
%      .orientation   - '111', '110', or '100' crystal plane normal
%                       aligned with the magnetic field
%
%      .nitrogen      - '14N' or '15N'
%
% Outputs:
%
%   sys   - Spinach system specification structure
%
%   inter - Spinach interaction specification structure
%
% alexey.bogdanov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=diamond_p1_13c.m>

function [sys,inter]=diamond_p1_13c(parameters)

% Check input count
if nargin~=1
    error('exactly one input argument is required.');
end

% Check consistency
grumble(parameters);

% Set crystal-to-defect frame rotation
crys2def=rotmat_align([1 1 1],[0 0 1]);

% Set the same orientation convention as diamond_p1.m
switch parameters.orientation

    case '100'

        R=rotmat_align([1 1 1],[1 0 0]);

    case '110'

        R=rotmat_align([1 1 1],[1 1 0]);

    case '111'

        R=eye(3);

    otherwise

        % Complain and bomb out
        error('unknown orientation specification.');

end

% Define the 13C site labels
site_lbl={'G1_C1','G2_C3','G3_calc','G4_C5','G5_calc',...
          'G6_calc','G7_calc','G8_C2','G9_calc','G10_calc',...
          'G11_calc','G12_calc','G13_calc','G14_C4',...
          'G15_calc','G16_calc','G17_calc','G18_calc'};

% Define spin system isotopes and labels
sys.isotopes=[{'E',parameters.nitrogen},repmat({'13C'},1,numel(site_lbl))];
sys.labels=[{'E',['N_' parameters.nitrogen]},site_lbl];
inter.zeeman.matrix=cell(1,numel(sys.isotopes));
inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes));

% Electron g-tensor
inter.zeeman.matrix{1}=R*diag([2.00220 2.00220 2.00218])*R';

% HFC and NQI tensor
switch parameters.nitrogen

    case '14N'

        % Hyperfine and quadrupolar coupling
        inter.coupling.matrix{1,2}=R*diag([81.3e6 81.3e6 114.0e6])*R';
        inter.coupling.matrix{2,2}=R*zfs2mat(-3.97e6,0,0,0,0)*R';

    case '15N'

        % Only hyperfine coupling, opposite sign
        inter.coupling.matrix{1,2}=R*diag([-114.0e6 -114.0e6 -159.9e6])*R';

    otherwise

        % Complain and bomb out
        error('wrong nitrogen isotope.');

end

% 13C site source map, in row order below:
%
% G1_C1    - C1, Cox et al. 1994, Table 2
% G2_C3    - C3, Cox et al. 1994, Table 2, with sign from Peaker et al. 2016
% G3_calc  - Peaker et al. 2016, Table 4
% G4_C5    - C5, Cox et al. 1994, Table 2, assignment from Peaker et al. 2016
% G5_calc  - Peaker et al. 2016, Table 4
% G6_calc  - Peaker et al. 2016, Table 4
% G7_calc  - Peaker et al. 2016, Table 4
% G8_C2    - C2, Cox et al. 1994, Table 2, assignment from Peaker et al. 2016
% G9_calc  - Peaker et al. 2016, Table 4
% G10_calc - Peaker et al. 2016, Table 4
% G11_calc - Peaker et al. 2016, Table 4
% G12_calc - Peaker et al. 2016, Table 4
% G13_calc - Peaker et al. 2016, Table 4
% G14_C4   - C4, Cox et al. 1994, Table 2, assignment from Peaker et al. 2016
% G15_calc - Peaker et al. 2016, Table 4
% G16_calc - Peaker et al. 2016, Table 4
% G17_calc - Peaker et al. 2016, Table 4
% G18_calc - Peaker et al. 2016, Table 4

% 13C hyperfine eigenvalues in MHz
hfc_vals=[+139.531 +139.531 +338.171;
          -26.488  -22.771  -25.319;
          -5.0     -4.4     -1.8;
          +11.757  +8.579   +8.122;
          +1.0     +1.1     +2.8;
          +3.7     +3.9     +5.8;
          +2.5     +2.7     +3.8;
          +30.921  +40.292  +31.662;
          +3.2     +3.5     +4.8;
          +2.7     +2.7     +4.8;
          +2.3     +2.4     +3.4;
          +1.0     +1.0     +1.7;
          +1.5     +1.6     +2.2;
          +10.638  +14.153  +10.618;
          +4.3     +4.3     +5.6;
          +0.8     +0.8     +1.4;
          +1.5     +1.5     +2.6;
          +1.2     +1.3     +1.8];

% 13C hyperfine principal-axis polar angles, degrees
hfc_theta=[90.00 35.264 54.736;
           90.00 52.36  37.64;
           90.00 13.00  77.00;
           71.50 138.60 55.00;
           54.00 37.00  85.00;
           90.00 59.00  31.00;
           90.00 75.00  15.00;
           90.00 58.66  31.34;
           90.00 23.00  67.00;
           71.00 90.00  19.00;
           90.00 68.00  22.00;
           90.00 52.00  38.00;
           90.00 64.00  26.00;
           90.00 59.19  30.81;
           90.00 29.00  61.00;
           35.00 90.00  55.00;
           90.00 29.00  61.00;
           28.00 90.00  62.00];

% 13C hyperfine principal-axis azimuths, degrees
hfc_phi=[135.00 225.00 45.00;
         315.00 45.00  225.00;
         315.00 45.00  225.00;
         33.20  101.00 137.00;
         21.00  211.00 115.00;
         315.00 45.00  225.00;
         135.00 45.00  225.00;
         315.00 45.00  225.00;
         135.00 225.00 45.00;
         224.00 315.00 45.00;
         315.00 225.00 45.00;
         315.00 45.00  225.00;
         135.00 45.00  225.00;
         315.00 45.00  225.00;
         135.00 45.00  225.00;
         225.00 315.00 45.00;
         135.00 225.00 45.00;
         45.00  135.00 225.00];

% Representative ideal-lattice directions for G1...G18
site_frac=[+0.25 +0.25 +0.25;
           +0.00 +0.50 +0.50;
           -0.25 -0.25 +0.25;
           +0.75 -0.25 +0.25;
           -0.50 +0.00 +0.50;
           +0.25 -0.75 +0.25;
           -0.75 +0.25 +0.25;
           +0.25 +0.75 +0.75;
           +0.00 +0.00 +1.00;
           +0.50 +0.50 +1.00;
           +0.25 +0.25 +1.25;
           +1.25 -0.25 -0.25;
           -0.25 -0.25 +1.25;
           +0.00 +1.00 +1.00;
           -1.00 +0.50 +0.50;
           +1.00 +1.00 +1.00;
           -1.00 -0.50 -0.50;
           -1.00 -1.00 +0.00];

% Peaker et al. tabulated distances from N-G1 bond mid-point, units of a0
site_rad=[0.284 0.545 0.553 0.735 0.735 0.894 0.894 0.899 0.906,...
          1.024 1.033 1.242 1.242 1.246 1.252 1.514 1.596 1.947];

% Set the diamond lattice parameter and the N-G1 bond mid-point
a0=3.567;
mid_frac=0.284*[1 1 1]/sqrt(3);

% Preallocate 13C coordinates
xyz_cub=zeros(numel(site_lbl),3);

% Reconstruct 13C coordinates from direction and radius
for n=1:numel(site_lbl)

    % Get radial direction from the N-G1 bond mid-point
    site_dir=site_frac(n,:)-mid_frac;
    site_dir=site_dir/norm(site_dir,2);

    % Put the site at the tabulated Peaker radius
    xyz_cub(n,:)=a0*(mid_frac+site_rad(n)*site_dir);

end

% Preallocate coordinate cell array
inter.coordinates=cell(1,numel(sys.isotopes));

% Leave the electron coordinate unknown
inter.coordinates{1}=[];

% Place nitrogen at the origin
inter.coordinates{2}=[0 0 0];

% Rotate coordinates into the requested orientation
xyz_def=(crys2def*xyz_cub')';
xyz_lab=(R*xyz_def')';
for n=1:numel(site_lbl)
    inter.coordinates{n+2}=xyz_lab(n,:);
end

% Build and rotate all 13C hyperfine tensors
for n=1:numel(site_lbl)

    % Convert principal-axis directions from polar angles
    axis_xyz=[sind(hfc_theta(n,:)).*cosd(hfc_phi(n,:));...
              sind(hfc_theta(n,:)).*sind(hfc_phi(n,:));...
              cosd(hfc_theta(n,:))];

    % Build the crystal-frame hyperfine tensor
    hfc_cub=axis_xyz*diag(1e6*hfc_vals(n,:))*axis_xyz';
    hfc_cub=(hfc_cub+hfc_cub')/2;

    % Rotate into the requested orientation
    hfc_def=crys2def*hfc_cub*crys2def';
    inter.coupling.matrix{1,n+2}=R*hfc_def*R';

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
if ~isfield(parameters,'nitrogen')
    error('parameters.nitrogen field is required.');
end
if(~ischar(parameters.nitrogen))
    error('parameters.nitrogen must be a character string.');
end
if ~ismember(parameters.nitrogen,{'14N','15N'})
    error('parameters.nitrogen must be ''14N'' or ''15N''.');
end
end

% Staff email at the University of Southampton, 09/12/2024:
% --------------------------------------------------------
%
%   Please be advised that an internal environmental audit [...] will take
%   place today at 10:00 am. The audit will be conducted by the following 
%   team members:
%
%       Dave Willard,    Senior Quality and Compliance Manager
%	    Sarah Puckett,   Environment and Sustainability Manager
%	    Wingki Tang,     Quality Management System Administrator
%	    Siobhan Balfour, Programme Administrator
% 
%   I [facilities manager] will be accompanying the team throughout all four 
%   buildings and assist with gathering the necessary data.
%
% Staff email at the University of Southampton, 21/11/2025:
% --------------------------------------------------------
%
%   Yesterday you will have received a communication from our Vice-Chancellor
%   outlining the challenges and pressures we are now facing as an institution
%   and the cost saving measures we need to apply, including a Voluntary Seve-
%   rance (VS) scheme. I appreciate that hearing this news will be unsettling
%   and for some will come as a surprise.

