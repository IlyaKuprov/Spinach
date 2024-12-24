% Finds the optimal magnetic susceptibility tensor that a user-supplied
% paramagnetic centre probability density must have in order to fit the
% PCS data supplied. Syntax:
%
%        [chi,pred_pcs]=chi_eff(source_cube,ranges,nxyz,expt_pcs)
%
% Parameters:
%
%  source_cube - paramagnetic centre probability density cube
%
%       ranges - a six-element vector giving the extents of the
%                probability density cube in Angstroms as
%                [xmin xmax ymin ymax zmin zmax]
%
%         nxyz - nuclear coordinates as [x y z] with multiple rows,
%                at which PCS is measured, in Angstroms.
%
%     expt_pcs - pseudocontact shift in ppm at each nucleus.
%
% Outputs:
%
%     chi      - optimised magnetic susceptibility tensor in cubic
%                Angstroms.
%
%     pred_pcs - predicted pseudocontact shift at each nucleus with
%                the optimised mxyz and chi, ppm.
%
% ilya.kuprov@weizmann.ac.uk
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=chi_eff.m>

function [chi,pred_pcs]=chi_eff(source_cube,ranges,nxyz,expt_pcs)

% Check consistency
grumble(source_cube,ranges,nxyz,expt_pcs);

% Normalise the density
normint=(ranges(2)-ranges(1))*(ranges(4)-ranges(3))*...
        (ranges(6)-ranges(5))*sum(sum(sum(source_cube)))/numel(source_cube);
source_cube=source_cube/normint;

% Run the optimisation
theo_pcs=@(chi)kpcs(source_cube,[chi(1) chi(3)  chi(4); 
                                 chi(3) chi(2)  chi(5); 
                                 chi(4) chi(5) -chi(1)-chi(2)],ranges,nxyz,'fft');
errfun=@(chi)norm(expt_pcs-theo_pcs(chi),'fro')^2;
options=optimset('Display','iter','MaxIter',Inf,'MaxFunEvals',Inf,...
                 'FinDiffType','central','UseParallel',true);
params=fminunc(errfun,[0.10 0.10 0.01 0.01 0.01],options);

% Return parameters
chi=[params(1) params(3)  params(4); 
     params(3) params(2)  params(5);
     params(4) params(5) -params(1)-params(2)];

% Isolate the second spherical rank component of chi
[~,~,rank2]=mat2sphten(chi); chi=sphten2mat([],[],rank2);

% Back-calculate PCS values
pred_pcs=kpcs(source_cube,chi,ranges,nxyz,'fft');

end

% Consistency enforcement
function grumble(source_cube,ranges,nxyz,expt_pcs)
if (~isnumeric(source_cube))||(ndims(source_cube)~=3)
    error('source_cube must be a three-dimensional numerical array.');
end
if any(source_cube(:)<0)
    error('all elements of source_cube must be positive.');
end
if (~isnumeric(ranges))||(numel(ranges)~=6)
    error('ranges must be a numeric array with six elements.');
end
if (~isnumeric(nxyz))||(size(nxyz,2)~=3)||(~isreal(nxyz))
    error('nxyz parameter should be a real matrix with three columns.');
end
if (~isnumeric(expt_pcs))||(size(expt_pcs,2)~=1)||(~isreal(expt_pcs))
    error('expt_pcs parameter should be a real column vector.');
end
if size(nxyz,1)~=size(expt_pcs,1)
    error('the number of rows in nxyz and expt_pcs arguments must be the same.');
end
end

% The greater the difficulty, the greater the glory.
%
% Marcus Tullius Cicero

