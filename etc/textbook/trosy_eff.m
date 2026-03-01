% TROSY efficiency in a two-spin system. Returns the extent
% of the cancellation of the CSA contribution to the trans-
% verse relaxation by the DD-CSA cross-correlation. Syntax:
%
%             eff=trosy_eff(B0,isotopes,xyz,csa)
%
% Parameters:
%
%    B0       - magnet field, Tesla
%
%    isotopes - a cell array with two character 
%               strings, e.g. {'19F','13C'} spe-
%               cifying spin-1/2 isotopes
%
%    xyz      - a cell array with two Cartesian
%               coordinate vectors in angstrom,
%               giving the locations of the two
%               nuclei
%
%    csa      - 3x3 chemical shielding or chemi-
%               cal shift (does not matter here)
%               tensor of the first spin in ppm;
%               its isotropic part will be drop-
%               ped automatically
%
% Outputs:
%
%    eff      - fraction of the CSA line width
%               that is compensated by DD-CSA
%               cross-correlation, 1 means that
%               the line width is purely dipolar
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=trosy_eff.m>

function eff=trosy_eff(B0,isotopes,xyz,csa)

% Check consistency
grumble(B0,isotopes,xyz,csa);

% Dipole-dipole coupling tensor
[~,~,~,~,DD]=xyz2dd(xyz{1},xyz{2},isotopes{1},isotopes{2});

% Anisotropic part of the Zeeman tensor at our field
Z=1e-6*B0*(csa-eye(3)*trace(csa)/3)*spin(isotopes{1});

% Second rank Blicharski invariants
[~,DsqZ]=blinv(Z); [~,X_DD_Z]=blprod(DD,Z);

% TROSY efficiency
eff=abs(X_DD_Z/DsqZ);

end

% Consistency enforcement
function grumble(B0,isotopes,xyz,csa)
if (~isnumeric(B0))||(~isreal(B0))||(~isscalar(B0))
    error('B0 must be a real number.');
end
if (~iscell(isotopes))||(numel(isotopes)~=2)||...
   (~ischar(isotopes{1}))||(~ischar(isotopes{2}))
    error('isotopes must be a cell array containing two character strings.');
end
if (~iscell(xyz))||(numel(xyz)~=2)||...
   (~isnumeric(xyz{1}))||(~isreal(xyz{1}))||(numel(xyz{1})~=3)||...
   (~isnumeric(xyz{2}))||(~isreal(xyz{2}))||(numel(xyz{2})~=3)
    error('xyz must be a cell array containing two three-element real vectors.');
end
if ~all(size(xyz{1})==size(xyz{2}))
    error('the two coordinate vectors in xyz must have the same dimension.');
end
if (~isnumeric(csa))||(~isreal(csa))||...
   (~ismatrix(csa))||any(size(csa)~=[3 3])
    error('csa must be a real 3x3 matrix.');
end
[~,mult_a]=spin(isotopes{1}); [~,mult_b]=spin(isotopes{2});
if (mult_a~=2)||(mult_b~=2)
    error('this function only supports spin-1/2 isotopes.');
end
if norm(xyz{1}-xyz{2},2)<0.25
    error('unphysical proximity.');
end
end

% Of caffeine's almighty healing powers 
% Which soothe away a thousand whisky sours,
% Transform the darkest morning into light,
% And bid the fiercest crapula take flight,
% I sing. Come Muse! Inspire my dreary tale
% That all may learn a lesson, none may fail.
%
% When Sol of morning brightness shone the ray
% Which pierced my curtains, heralding new day,
% I rose up and, stumbling from the bed,
% Perceived a battle raging in my head
% Between some Scottish malt and Spanish cava:
% I knew at once the only cure was Java.
%
% The kitchen soon I gained on falt'ring feet,
% Gazed on the jars of storage - shining, neat -
% And saw forthwith my life was quite undone.
% Was there Nespresso? Reader, there was none.
%
% Tom Adam

