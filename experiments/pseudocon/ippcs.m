% Fits the point electron model PCS to the experimental pseudocon-
% tact shift coordinates and values. Syntax:
%
%       [exyz,chi,pred_pcs]=ippcs(nxyz,mguess,expt_pcs)
%
% Parameters:
%
%     nxyz     - nuclear coordinates as [x y z] with multiple rows,
%                at which PCS is measured, in Angstroms.
%
%     mguess   - initial guess for the unpaired electron coordina-
%                tes as [x y z], in Angstroms.
%
%     expt_pcs - pseudocontact shift in ppm at each nucleus.
%
% Outputs:
%
%     mxyz     - optimized paramagnetic centre coordinates as [x y z],
%                in Angstroms.
%
%     chi      - optimized magnetic susceptibility tensor in cubic
%                Angstroms.
%
%     pred_pcs - predicted pseudocontact shift at each nucleus with
%                the optimized mxyz and chi, ppm.
%
%     s_mxyz   - standard deviations of paramagnetic centre 
%                coordinates as [x y z], in Angstroms.
%
%     s_chi    - standard deviations of magnetic susceptibility 
%                tensor elements in cubic Angstroms.
%
% Note: a good initial guess for the paramagnetic centre location is
%       essential for a successful fit.
% 
% ilya.kuprov@weizmann.ac.il
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=Ippcs.m>

function [mxyz,chi,pred_pcs,s_mxyz,s_chi]=ippcs(nxyz,mguess,expt_pcs)

% Check consistency
grumble(nxyz,mguess,expt_pcs);

% Set optimisation options
warning('off','optim:fminunc:SwitchingMethod');
options=optimset('Display','iter','MaxIter',inf,'UseParallel',true,...
                 'MaxFunEvals',inf,'FinDiffType','central');

% Vector of squared residuals
vec_res_sq=@(x)(expt_pcs-ppcs(nxyz,x(1:3),x(4:8))).^2;

% Sum of squared residuals 
sum_res_sq=@(x)sum(vec_res_sq(x));

% Fit the point model
p=fminunc(sum_res_sq,[mguess 0.1 0.1 0.1 0.1 0.1],options);

% Get the predicted PCS values
pred_pcs=ppcs(nxyz,p(1:3),p(4:8));

% Compute the Jacobian at the optimal point
jac=jacobianest(vec_res_sq,p);

% Get the Studentized residual
sdr=sqrt(sum_res_sq(p)/(numel(expt_pcs)-8));

% Get the standard deviations
sp=sqrt(diag((sdr^2)*inv(jac'*jac)))'; %#ok<MINV>

% Return parameters and standard deviations
mxyz=p(1:3); s_mxyz=sp(1:3); 
chi=[p(4)     p(5)     p(6); 
     p(5)     p(7)     p(8); 
     p(6)     p(8) -p(4)-p(7)];
s_chi=[sp(4)          sp(5)          sp(6);
       sp(5)          sp(7)          sp(8);
       sp(6)          sp(8)  sqrt(sp(4)^2+sp(7)^2)];

end

% Consistency enforcement
function grumble(nxyz,mguess,expt_pcs)
if iscell(nxyz)
    for n=1:numel(nxyz)
        if (~isnumeric(nxyz{n}))||(size(nxyz{n},2)~=3)||(~isreal(nxyz{n}))
            error('nxyz parameter should be a cell array of real matrices with three columns.');
        end
        if size(nxyz{n},1)~=size(expt_pcs,1)
            error('the number of rows in nxyz and expt_pcs arguments must be the same.');
        end
    end
else
    if (~isnumeric(nxyz))||(size(nxyz,2)~=3)||(~isreal(nxyz))
        error('nxyz parameter should be a real matrix with three columns.');
    end
    if size(nxyz,1)~=size(expt_pcs,1)
        error('the number of rows in nxyz and expt_pcs arguments must be the same.');
    end
end
if (~isnumeric(mguess))||(size(mguess,2)~=3)||(size(mguess,1)~=1)||(~isreal(mguess))
    error('mguess parameter should be a real row vector with three elements.');
end
if (~isnumeric(expt_pcs))||(size(expt_pcs,2)~=1)||(~isreal(expt_pcs))
    error('expt_pcs parameter should be a real column vector.');
end
end

% The more one is hated, I find, the happier one is. 
%
% Louis Ferdinand Celine

