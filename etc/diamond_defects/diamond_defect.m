% Spin system for EPR-active defects in diamond. Syntax:
%
%          [sys,inter]=diamond_defect(parameters)
%
% This is a literature-parameter database for diamond defects
% other than the standard NV- ground-state and P1 entries.
%
% Parameters:
%
%    parameters.defect      - defect label: 'nv0_es', 'n2vm',
%                             'war9', 'war10', 'r2', 'r4_w6',
%                             'w29', 'r5', 'o1', 'r6', 'r10',
%                             'r11', 'siv0',
%                             'gev0', 'w8', 'ne1', 'ne2',
%                             'ne3', 'ne4', 'ne5', 'ne8',
%                             'ab1', 'ab2', 'ab3', 'ab4',
%                             'ab5', 'nol1', 'n3', 'ok1',
%                             'o4', 'nlo2', 'ma1', 'np1',
%                             'np2', 'np3', 'np4', 'np5',
%                             'np6', 'np8', or 'np9'
%
%    parameters.orientation - crystal plane normal aligned with
%                             the magnetic field: '111', '110',
%                             or '100', default is '111'
%
%    parameters.nitrogen    - '14N' or '15N' where supported,
%                             default is defect-specific
%
%    parameters.include_13c - include reported 13C hyperfine
%                             couplings where available, false
%                             by default
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% Sources:
%
%    Felton et al., Phys. Rev. B 77, 081201 (2008),
%    https://doi.org/10.1103/PhysRevB.77.081201
%    Green et al., Phys. Rev. B 92, 165204 (2015),
%    https://doi.org/10.1103/PhysRevB.92.165204
%    Felton et al., J. Phys. Condens. Matter 21, 364212 (2009),
%    https://doi.org/10.1088/0953-8984/21/36/364212
%    Hunt et al., Phys. Rev. B 61, 3863 (2000),
%    https://doi.org/10.1103/PhysRevB.61.3863
%    Kirui et al., Diam. Relat. Mater. 8, 1569 (1999),
%    https://doi.org/10.1016/S0925-9635(99)00037-0
%    Iakoubovskii and Stesmans, Phys. Rev. B 66, 045406 (2002),
%    https://doi.org/10.1103/PhysRevB.66.045406
%    Edmonds et al., Phys. Rev. B 77, 245205 (2008),
%    https://doi.org/10.1103/PhysRevB.77.245205
%    Nadolinny et al., Phys. Status Solidi A 213, 2623 (2016),
%    https://doi.org/10.1002/pssa.201600211
%    Ludwig and Woodbury, Phys. Rev. B 41, 3905 (1990),
%    https://doi.org/10.1103/PhysRevB.41.3905
%    Ball, PhD thesis, OIST Graduate University (2021)
%    Nadolinny et al., Crystals 7, 237 (2017),
%    https://doi.org/10.3390/cryst7080237
%
% Notes:
%
%    Hyperfine and zero-field splitting constants reported in
%    magnetic-field units are converted using the free-electron
%    gyromagnetic ratio. Where the source gives unresolved signs
%    for D, the positive sign is used unless parameters.d_sign is
%    supplied.
%
% <https://spindynamics.org/wiki/index.php?title=diamond_defect.m>

function [sys,inter]=diamond_defect(parameters)

% Check consistency
grumble(parameters);

% Set default parameters
if ~isfield(parameters,'orientation')
    parameters.orientation='111';
end
if ~isfield(parameters,'include_13c')
    parameters.include_13c=false;
end
if ~isfield(parameters,'d_sign')
    parameters.d_sign=1;
end

% Set useful constants
hz_per_mt=abs(spin('E'))/(2*pi)*1e-3;
hz_per_t=1e3*hz_per_mt;

% Set crystal-to-laboratory rotation
C=cryst_rot(parameters.orientation);

% Initialise the spin specification
nuclei={}; zfs=[];
defect=lower(parameters.defect);

% Select the defect specification
switch defect

    case {'nv0_es','nv0'}

        % Excited quartet state of the neutral NV centre
        % DOI: https://doi.org/10.1103/PhysRevB.77.081201
        electron='E4'; frame=frame_z([1 1 1]);
        gmat=tensor([2.0035 2.0035 2.0029],frame);
        zfs=tensor_zfs(1685e6,0,frame);
        nitrogen=getfield_def(parameters,'nitrogen','15N');
        if ~strcmp(nitrogen,'15N')
            error('NV0 excited-state parameters are only available here for 15N.');
        end
        nuclei=add_nuc(nuclei,'15N',tensor([-23.8e6 -23.8e6 -35.7e6],frame),[]);

    case {'n2vm','n2v-'}

        % Negatively charged N2V centre
        % DOI: https://doi.org/10.1103/PhysRevB.92.165204
        electron='E'; c2rot=rot_axis([0 0 1],180);
        gframe=frame_xyz([1 1 0],[0 0 1],[1 -1 0]);
        nframe=rot_axis([1 -1 0],-3.5)*...
               frame_xyz([1 1 -2],[1 1 1],[1 -1 0]);
        cframe=rot_axis([-1 -1 0],2.0)*...
               frame_xz([-1 -1 0],[-1 1 1]);
        gmat=tensor([2.00345 2.00274 2.00271],gframe);
        nitrogen=getfield_def(parameters,'nitrogen','15N');
        if strcmp(nitrogen,'15N')
            nuclei=add_nuc(nuclei,'15N',tensor([3.47e6 4.51e6 4.09e6],nframe),[]);
            nuclei=add_nuc(nuclei,'15N',tensor([3.47e6 4.51e6 4.09e6],c2rot*nframe),[]);
        elseif strcmp(nitrogen,'14N')
            scale=spin('14N')/spin('15N');
            Qframe=frame_z([1 1 1]);
            nuclei=add_nuc(nuclei,'14N',scale*tensor([3.47e6 4.51e6 4.09e6],nframe),...
                           tensor_zfs(-5.0e6,0,Qframe));
            nuclei=add_nuc(nuclei,'14N',scale*tensor([3.47e6 4.51e6 4.09e6],c2rot*nframe),...
                           tensor_zfs(-5.0e6,0,c2rot*Qframe));
        else
            error('parameters.nitrogen must be ''14N'' or ''15N''.');
        end
        if parameters.include_13c
            nuclei=add_nuc(nuclei,'13C',tensor([202.3e6 202.3e6 317.5e6],cframe),[]);
            nuclei=add_nuc(nuclei,'13C',tensor([202.3e6 202.3e6 317.5e6],c2rot*cframe),[]);
        end

    case {'war9','war10'}

        % Nitrogen split-interstitial defects
        % DOI: https://doi.org/10.1088/0953-8984/21/36/364212
        electron='E'; nitrogen=getfield_def(parameters,'nitrogen','15N');
        if ~strcmp(nitrogen,'15N')
            error('WAR9/WAR10 parameters are only available here for 15N.');
        end
        gframe=frame_xyz(sph_vec(90,45),sph_vec(180,45),sph_vec(90,315));
        if strcmp(defect,'war9')
            gmat=tensor([2.00343 2.00272 2.00268],gframe);
            Amat=tensor([8.30e6 7.85e6 8.17e6],gframe);
        else
            aframe=frame_xyz(sph_vec(44.8,45.0),sph_vec(134.8,45.0),sph_vec(90,315));
            gmat=tensor([2.00344 2.00272 2.00269],gframe);
            Amat=tensor([1.00e6 -1.01e6 0.00e6],aframe);
        end
        nuclei=add_nuc(nuclei,'15N',Amat,[]);

    case 'r2'

        % Neutral <100>-split self-interstitial
        % DOI: https://doi.org/10.1103/PhysRevB.61.3863
        electron='E3'; frame=frame_z([1 0 0]);
        gmat=tensor([2.0019 2.0019 2.0021],frame);
        zfs=tensor_zfs(parameters.d_sign*4173e6,0,frame);

    case {'r4_w6','w6','r4','w29','r5','o1','r6','r10','r11'}

        % Vacancy, divacancy, and vacancy-chain centres
        [electron,giso,zfs]=vacancy_data(defect);
        gmat=eye(3)*giso;

    case 'siv0'

        % Neutral silicon-vacancy centre
        % DOI: https://doi.org/10.1103/PhysRevB.77.245205
        electron='E3'; frame=frame_z([1 1 1]);
        gmat=tensor([2.0035 2.0035 2.0042],frame);
        zfs=tensor_zfs(1000e6,0,frame);
        silicon=getfield_def(parameters,'silicon','29Si');
        if strcmp(silicon,'29Si')
            nuclei=add_nuc(nuclei,'29Si',tensor([78.9e6 78.9e6 76.3e6],frame),[]);
        elseif ~strcmp(silicon,'none')
            nuclei=add_nuc(nuclei,silicon,zeros(3),[]);
        end
        if parameters.include_13c
            Cmat=tensor([30.2e6 30.2e6 66.2e6],frame);
            for n=1:6
                nuclei=add_nuc(nuclei,'13C',Cmat,[]);
            end
        end

    case 'gev0'

        % Neutral germanium-vacancy centre
        % DOI: https://doi.org/10.1002/pssa.201600211
        electron='E3'; frame=frame_z([1 1 1]);
        gmat=tensor([2.0027 2.0027 2.0025],frame);
        zfs=tensor_zfs(80.3*hz_per_mt,0,frame);
        germanium=getfield_def(parameters,'germanium','73Ge');
        if strcmp(germanium,'73Ge')
            nuclei=add_nuc(nuclei,'73Ge',eye(3)*1.64*hz_per_mt,[]);
        elseif ~strcmp(germanium,'none')
            nuclei=add_nuc(nuclei,germanium,zeros(3),[]);
        end

    case 'w8'

        % Substitutional nickel W8 centre
        % DOI: https://doi.org/10.1103/PhysRevB.41.3905
        electron='E4';
        gmat=eye(3)*2.032;
        nickel=getfield_def(parameters,'nickel','61Ni');
        if strcmp(nickel,'61Ni')
            nuclei=add_nuc(nuclei,'61Ni',eye(3)*0.65*hz_per_mt,[]);
        elseif ~strcmp(nickel,'none')
            nuclei=add_nuc(nuclei,nickel,zeros(3),[]);
        end
        if parameters.include_13c
            Cmat=tensor([0.340 0.340 1.339]*hz_per_mt,frame_z([1 1 1]));
            for n=1:4
                nuclei=add_nuc(nuclei,'13C',Cmat,[]);
            end
        end

    case {'ne1','ne2','ne3','ne5','ne8'}

        % Nickel-nitrogen split-vacancy centres
        % DOI: https://doi.org/10.3390/cryst7080237
        electron='E'; [gvals,avalues,alpha]=nickel_ne_data(defect);
        frame=frame_alpha(alpha);
        gmat=tensor(gvals,frame);
        for n=1:size(avalues,1)
            nuclei=add_nuc(nuclei,'14N',tensor(avalues(n,:)*hz_per_mt,frame),[]);
        end

    case 'ne4'

        % Nickel split-vacancy NE4 centre
        % DOI: https://doi.org/10.3390/cryst7080237
        electron='E'; frame=frame_z([1 1 1]);
        gmat=tensor([2.0988 2.0988 2.0227],frame);

    case {'ab1','ab2','ab3','ab4'}

        % Nitrogen-free nickel centres without reported HFS
        % DOI: https://doi.org/10.3390/cryst7080237
        electron='E'; [gvals,frame]=nickel_ab_data(defect);
        gmat=tensor(gvals,frame);

    case 'ab5'

        % Triplet nitrogen-free nickel AB5 centre
        % DOI: https://doi.org/10.3390/cryst7080237
        electron='E3'; frame=frame_z([1 1 1]);
        gmat=tensor([2.022 2.022 2.037],frame);
        zfs=tensor_zfs(parameters.d_sign*1.132*hz_per_t,0,frame);

    case {'nol1','nirim5'}

        % Triplet nitrogen-free nickel NOL1/NIRIM5 centre
        % DOI: https://doi.org/10.3390/cryst7080237
        electron='E3'; frame=frame_z([1 1 1]);
        gmat=tensor([2.002 2.002 2.0235],frame);
        zfs=tensor_zfs(-6.10*hz_per_t,0,frame);

    case {'n3','ok1'}

        % Titanium-containing centres
        % DOI: https://doi.org/10.3390/cryst7080237
        electron='E'; [gvals,An,Ati,g_alpha,A_alpha]=titanium_data(defect);
        gframe=frame_alpha(g_alpha); aframe=frame_alpha(A_alpha);
        gmat=tensor(gvals,gframe);
        nuclei=add_nuc(nuclei,'14N',tensor(An*hz_per_mt,aframe),[]);
        titanium=getfield_def(parameters,'titanium','47Ti');
        if ~strcmp(titanium,'none')
            nuclei=add_nuc(nuclei,titanium,tensor(Ati*hz_per_mt,aframe),[]);
        end
        if strcmp(defect,'ok1')&&parameters.include_13c
            Cmat=tensor([2.62 2.62 4.38]*hz_per_mt,frame_xz([1 1 0],[1 -1 -1]));
            nuclei=add_nuc(nuclei,'13C',Cmat,[]);
            Cmat=tensor([2.62 2.62 4.38]*hz_per_mt,frame_xz([1 1 0],[-1 1 -1]));
            nuclei=add_nuc(nuclei,'13C',Cmat,[]);
        end

    case {'o4','nlo2'}

        % Cobalt-containing centres
        % DOI: https://doi.org/10.3390/cryst7080237
        electron='E'; [gvals,Aco,alpha]=cobalt_data(defect);
        frame=frame_cobalt(alpha);
        gmat=tensor(gvals,frame);
        nuclei=add_nuc(nuclei,'59Co',tensor(Aco*hz_per_mt,frame),[]);

    case {'ma1','np1','np2','np3','np4','np5','np6','np8','np9'}

        % Phosphorus-containing centres
        % DOI: https://doi.org/10.3390/cryst7080237
        electron='E'; [gmat,nuclei]=phosphorus_data(defect,hz_per_mt,parameters.include_13c);

    otherwise

        % Complain and bomb out
        error('unknown diamond defect specification.');

end

% Build the Spinach structures
sys.isotopes={electron};
inter.zeeman.matrix{1}=C*gmat*C';
if ~isempty(zfs)
    inter.coupling.matrix{1,1}=traceless(C*zfs*C');
else
    inter.coupling.matrix{1,1}=[];
end
for n=1:numel(nuclei)
    sys.isotopes{n+1}=nuclei{n}.iso;
    inter.zeeman.matrix{n+1}=zeros(3);
    inter.coupling.matrix{1,n+1}=C*nuclei{n}.A*C';
    if ~isempty(nuclei{n}.Q)
        inter.coupling.matrix{n+1,n+1}=traceless(C*nuclei{n}.Q*C');
    else
        inter.coupling.matrix{n+1,n+1}=[];
    end
end

end

% Field getter with a default value
function value=getfield_def(parameters,name,default)
if isfield(parameters,name)
    value=parameters.(name);
else
    value=default;
end
end

% Add a nucleus to the local specification
function nuclei=add_nuc(nuclei,iso,A,Q)
nucleus=struct('iso',iso,'A',A,'Q',Q);
nuclei{end+1}=nucleus;
end

% Build a tensor from principal values and axes
function M=tensor(values,frame)
frame=orth_frame(frame);
M=frame*diag(values)*frame';
M=(M+M')/2;
end

% Build a zero-field splitting tensor from principal axes
function M=tensor_zfs(D,E,frame)
frame=orth_frame(frame);
M=traceless(frame*zfs2mat(D,E,0,0,0)*frame');
end

% Enforce symmetric traceless form for quadratic couplings
function M=traceless(M)
M=(M+M')/2;
M=M-eye(3)*trace(M)/3;
M=(M+M')/2;
end

% Crystal-to-laboratory rotation matrix
function C=cryst_rot(orientation)
switch orientation
    case '111'
        C=rotmat_align([1 1 1],[0 0 1]);
    case '110'
        C=rotmat_align([1 1 0],[0 0 1]);
    case '100'
        C=rotmat_align([1 0 0],[0 0 1]);
    otherwise
        error('unknown orientation specification.');
end
end

% Make a principal-axis frame from the z axis
function frame=frame_z(zaxis)
zaxis=zaxis(:)/norm(zaxis,2);
if abs(dot(zaxis,[0;0;1]))<0.9
    xaxis=cross([0;0;1],zaxis);
else
    xaxis=cross([0;1;0],zaxis);
end
xaxis=xaxis/norm(xaxis,2);
yaxis=cross(zaxis,xaxis);
frame=[xaxis yaxis zaxis];
end

% Make a principal-axis frame from three vectors
function frame=frame_xyz(xaxis,yaxis,zaxis)
frame=[xaxis(:) yaxis(:) zaxis(:)];
frame=orth_frame(frame);
end

% Orthogonalise a right-handed frame
function frame=orth_frame(frame)
[frame,~]=qr(frame,0);
if det(frame)<0
    frame(:,3)=-frame(:,3);
end
end

% Vector from polar angles in the diamond crystal frame
function v=sph_vec(theta,phi)
v=[sind(theta)*cosd(phi);sind(theta)*sind(phi);cosd(theta)];
end

% Rotation matrix for an axis-angle rotation in degrees
function R=rot_axis(axis,angle)
axis=axis(:)/norm(axis,2);
K=[0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
R=eye(3)+sind(angle)*K+(1-cosd(angle))*(K*K);
end

% Principal-axis frame from a specified x axis and z axis
function frame=frame_xz(xaxis,zaxis)
xaxis=xaxis(:)/norm(xaxis,2);
zaxis=zaxis(:)-xaxis*dot(xaxis,zaxis(:));
zaxis=zaxis/norm(zaxis,2);
yaxis=cross(zaxis,xaxis);
frame=frame_xyz(xaxis,yaxis,zaxis);
end

% Frame used in the Nadolinny nickel and titanium tables
function frame=frame_alpha(alpha)
xaxis=[1;-1;0]/sqrt(2);
ybase=[1;1;0]/sqrt(2);
zbase=[0;0;1];
yaxis=cosd(alpha)*ybase+sind(alpha)*zbase;
zaxis=cross(xaxis,yaxis);
frame=frame_xyz(xaxis,yaxis,zaxis);
end

% Frame used in the Nadolinny cobalt table
function frame=frame_cobalt(alpha)
yaxis=[0;1;1]/sqrt(2);
xbase=[1;0;0];
zbase=cross(xbase,yaxis);
xaxis=cosd(alpha)*xbase+sind(alpha)*zbase;
zaxis=cross(xaxis,yaxis);
frame=frame_xyz(xaxis,yaxis,zaxis);
end

% Nickel-nitrogen centre table data
function [gvals,avalues,alpha]=nickel_ne_data(defect)
switch defect
    case 'ne1'
        % DOI: https://doi.org/10.3390/cryst7080237
        gvals=[2.1282 2.0070 2.0908]; alpha=14;
        avalues=[2.09 1.43 1.45; 2.09 1.43 1.45];
    case 'ne2'
        % DOI: https://doi.org/10.3390/cryst7080237
        gvals=[2.1301 2.0100 2.0931]; alpha=14;
        avalues=[2.10 1.42 1.41; 1.87 1.18 1.25; 0.18 0.35 0.25];
    case 'ne3'
        % DOI: https://doi.org/10.3390/cryst7080237
        gvals=[2.0729 2.0100 2.0476]; alpha=14;
        avalues=[1.60 1.24 1.15; 0.66 0.50 0.50; 0.66 0.50 0.50];
    case 'ne5'
        % DOI: https://doi.org/10.3390/cryst7080237
        gvals=[2.0329 2.0898 2.0476]; alpha=27.5;
        avalues=[1.22 0.98 0.89; 1.22 0.98 0.89];
    case 'ne8'
        % DOI: https://doi.org/10.3390/cryst7080237
        gvals=[2.0439 2.1722 2.0476]; alpha=27.5;
        avalues=[1.14 0.78 0.75; 1.14 0.78 0.75; 1.14 0.78 0.75; 1.14 0.78 0.75];
    otherwise
        error('unknown nickel NE centre.');
end
end

% Nitrogen-free nickel centre table data
function [gvals,frame]=nickel_ab_data(defect)
switch defect
    case 'ab1'
        % DOI: https://doi.org/10.3390/cryst7080237
        frame=frame_z([1 1 1]); gvals=[2.0920 2.0920 2.0024];
    case 'ab2'
        % DOI: https://doi.org/10.3390/cryst7080237
        frame=frame_z([1 1 1]); gvals=[2.0672 2.0672 2.0072];
    case 'ab3'
        % DOI: https://doi.org/10.3390/cryst7080237
        frame=frame_xyz([1 0 0],[0 1 1],[0 -1 1]); gvals=[2.1105 2.0663 2.0181];
    case 'ab4'
        % DOI: https://doi.org/10.3390/cryst7080237
        frame=frame_xyz([1 0 0],[0 1 1],[0 -1 1]); gvals=[2.0220 2.0094 2.0084];
    otherwise
        error('unknown nickel AB centre.');
end
end

% Vacancy and vacancy-chain centre table data
function [electron,giso,zfs]=vacancy_data(defect)
% The cited vacancy tables fix the ZFS tensors; the g values are
% isotropic near-free-electron values used for the corresponding
% EPR assignments, not resolved anisotropic g tensors.
switch defect
    case {'r4_w6','w6','r4'}
        % Neutral divacancy R4/W6, S=1, principal z approximately [111]
        electron='E3'; giso=2.0022;
        % DOI: https://doi.org/10.1016/S0925-9635(99)00037-0
        zfs=tensor([105 197 -303]*1e6,frame_z([1 1 1]));
    case 'w29'
        % Negative divacancy W29, S=3/2, principal z approximately [111]
        electron='E4'; giso=2.0019;
        % DOI: https://doi.org/10.1016/S0925-9635(99)00037-0
        zfs=tensor([297 156 -453]*1e6,frame_z([1 1 1]));
    case 'r5'
        % Three-vacancy chain R5, S=1, principal z along [110]
        electron='E3'; giso=2.0023;
        % DOI: https://doi.org/10.1103/PhysRevB.66.045406
        zfs=tensor([283 244 -524]*1e6,frame_z([1 1 0]));
    case 'o1'
        % Four-vacancy chain O1, S=1, principal z along [110]
        electron='E3'; giso=2.0023;
        % DOI: https://doi.org/10.1103/PhysRevB.66.045406
        zfs=tensor([109 95 -205]*1e6,frame_z([1 1 0]));
    case 'r6'
        % Five-vacancy chain R6, S=1, principal z along [110]
        electron='E3'; giso=2.0023;
        % DOI: https://doi.org/10.1103/PhysRevB.66.045406
        zfs=tensor([62 59 -120]*1e6,frame_z([1 1 0]));
    case 'r10'
        % Six-vacancy chain R10, S=1, principal z along [110]
        electron='E3'; giso=2.0023;
        % DOI: https://doi.org/10.1103/PhysRevB.66.045406
        zfs=tensor([36 36 -73]*1e6,frame_z([1 1 0]));
    case 'r11'
        % Seven-vacancy chain R11, S=1, principal z along [110]
        electron='E3'; giso=2.0023;
        % DOI: https://doi.org/10.1103/PhysRevB.66.045406
        zfs=tensor([27 27 -53]*1e6,frame_z([1 1 0]));
    otherwise
        error('unknown vacancy centre.');
end
end

% Titanium centre table data
function [gvals,An,Ati,g_alpha,A_alpha]=titanium_data(defect)
switch defect
    case 'n3'
        % DOI: https://doi.org/10.3390/cryst7080237
        gvals=[2.0022 2.0025 2.0020];
        An=[0.11 0.15 0.11];
        Ati=[0.28 0.40 0.28];
        g_alpha=32; A_alpha=26;
    case 'ok1'
        % DOI: https://doi.org/10.3390/cryst7080237
        gvals=[2.0031 2.0019 2.0025];
        An=[0.55 0.77 0.54];
        Ati=[0.06 0.06 0.06];
        g_alpha=40; A_alpha=20;
    otherwise
        error('unknown titanium centre.');
end
end

% Cobalt centre table data
function [gvals,Aco,alpha]=cobalt_data(defect)
switch defect
    case 'o4'
        % DOI: https://doi.org/10.3390/cryst7080237
        gvals=[2.3463 1.8438 1.7045];
        Aco=[8.86 6.43 5.82];
        alpha=29;
    case 'nlo2'
        % DOI: https://doi.org/10.3390/cryst7080237
        gvals=[2.3277 1.7982 1.7149];
        Aco=[8.24 6.57 5.76];
        alpha=28;
    otherwise
        error('unknown cobalt centre.');
end
end

% Phosphorus centre table data
function [gmat,nuclei]=phosphorus_data(defect,hz_per_mt,include_13c)
nuclei={}; frame=eye(3);
switch defect
    case 'ma1'
        % DOI: https://doi.org/10.3390/cryst7080237
        gmat=eye(3)*2.0025;
        nuclei=add_nuc(nuclei,'31P',tensor([1.96 1.96 2.32]*hz_per_mt,frame_z([1 1 1])),[]);
        if include_13c
            nuclei=add_nuc(nuclei,'13C',tensor([13.92 13.92 18.13]*hz_per_mt,frame_z([1 1 1])),[]);
        end
    case 'np1'
        % DOI: https://doi.org/10.3390/cryst7080237
        gmat=tensor([2.00243 2.0028 2.0026],frame);
        nuclei=add_nuc(nuclei,'31P',tensor([2.08 2.02 2.18]*hz_per_mt,frame),[]);
        nuclei=add_nuc(nuclei,'14N',tensor([4.08 3.10 3.00]*hz_per_mt,frame),[]);
    case 'np2'
        % DOI: https://doi.org/10.3390/cryst7080237
        gmat=eye(3)*2.0025;
        nuclei=add_nuc(nuclei,'31P',tensor([2.09 2.09 2.34]*hz_per_mt,frame_z([1 1 1])),[]);
        nuclei=add_nuc(nuclei,'14N',tensor([3.09 3.09 6.42]*hz_per_mt,frame_z([1 1 1])),[]);
    case 'np3'
        % DOI: https://doi.org/10.3390/cryst7080237
        gmat=eye(3)*2.0025;
        nuclei=add_nuc(nuclei,'31P',tensor([18.23 18.23 17.48]*hz_per_mt,frame_z([1 1 1])),[]);
        nuclei=add_nuc(nuclei,'14N',tensor([0.33 0.33 0.10]*hz_per_mt,frame_z([1 1 1])),[]);
    case 'np4'
        % DOI: https://doi.org/10.3390/cryst7080237
        gmat=tensor([2.0009 2.0012 2.00047],frame);
        nuclei=add_nuc(nuclei,'31P',tensor([5.456 3.838 3.80]*hz_per_mt,frame),[]);
    case 'np5'
        % DOI: https://doi.org/10.3390/cryst7080237
        gmat=tensor([2.0009 2.0009 2.00087],frame_z([1 1 1]));
        nuclei=add_nuc(nuclei,'31P',tensor([1.024 1.024 6.522]*hz_per_mt,frame_z([1 1 1])),[]);
    case 'np6'
        % DOI: https://doi.org/10.3390/cryst7080237
        gmat=tensor([2.00083 2.00083 2.00085],frame_z([1 1 1]));
        nuclei=add_nuc(nuclei,'31P',tensor([7.585 2.942 2.328]*hz_per_mt,frame),[]);
    case 'np8'
        % DOI: https://doi.org/10.3390/cryst7080237
        gmat=tensor([2.0016 2.0016 2.0048],frame_z([1 1 1]));
        nuclei=add_nuc(nuclei,'31P',tensor([3.2 3.2 5.6]*hz_per_mt,frame_z([1 1 1])),[]);
        nuclei=add_nuc(nuclei,'31P',tensor([8.8 8.8 13.6]*hz_per_mt,frame_z([1 1 1])),[]);
    case 'np9'
        % DOI: https://doi.org/10.3390/cryst7080237
        gmat=tensor([2.0038 2.0038 2.0030],frame_z([1 1 1]));
        nuclei=add_nuc(nuclei,'31P',tensor([2.2 2.2 1.4]*hz_per_mt,frame_z([1 1 1])),[]);
    otherwise
        error('unknown phosphorus centre.');
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
if isfield(parameters,'orientation')&&(~ischar(parameters.orientation))
    error('parameters.orientation must be a character string.');
end
if isfield(parameters,'nitrogen')&&(~ischar(parameters.nitrogen))
    error('parameters.nitrogen must be a character string.');
end
if isfield(parameters,'include_13c')&&(~islogical(parameters.include_13c))
    error('parameters.include_13c must be logical.');
end
if isfield(parameters,'d_sign')&&(~isnumeric(parameters.d_sign)||~isreal(parameters.d_sign)||~isscalar(parameters.d_sign))
    error('parameters.d_sign must be a real scalar.');
end
end

% A fact should be simple enough to survive repetition, but not
% so simple that repetition changes its meaning.

