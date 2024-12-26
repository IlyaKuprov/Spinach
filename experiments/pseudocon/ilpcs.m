% Fits experimental PCS data using the distributed paramagnetc centre 
% model described in
%
%                  http://dx.doi.org/10.1039/c6cp05437d
%
% Syntax:
%
%       [mxyz,chi,Ilm,pred_pcs,s_mxyz,s_chi,s_Ilm]=...
%                             ilpcs(nxyz,expt_pcs,ranks,mguess)
%
% Parameters: 
%
%      nxyz  - nuclear coordinates as [x y z] with multiple rows,
%              at which PCS is to be evaluated, in Angstroms.
%
%  expt_pcs  - a column vector of experimental pseudocontact shifts
%              in ppm  
%
%      ranks - row of multipole expansion ranks to be used in the
%              fitting procedure
%
%     mguess - guess value for the paramagnetic centre position,
%              a three-element vector in Angstrom
%             
% Output:
%
%     mxyz   - optimized paramagnetic centre coordinates as [x y z],
%              in Angstroms.
%
%     chi    - optimized magnetic susceptibility tensor in cubic
%              Angstroms.
%
%     Ilm    - {[],[]} cell array of numbers corresponding to the 
%              multipole moments defined in the paper cited above:
%              
%                for L=0,  Ilm=N/2/sqrt(pi)
%              
%                for L=1,  Ilm=[real(I11) I10 imag(I11)]
%              
%                for L=2,  Ilm=[real(I22) real(I21) I20 imag(I21) imag(I22)]
%
%              et cetera.
%
%  pred_pcs  - predicted pseudocontact shift (in ppm) at each of 
%              the nuclei.
%
%     chi    - optimized magnetic susceptibility tensor in cubic
%              Angstroms.
%
%     s_mxyz - standard deviations of paramagnetic centre 
%              coordinates as [x y z], in Angstroms.
%
%     s_chi  - standard deviations of magnetic susceptibility 
%              tensor elements in cubic Angstroms.
%
%     s_Ilm  - standard deviations of the multipole moments, arranged
%              in the same order as the moments themselves.
%
% Note: a good initial guess for the paramagnetic centre location is
%       essential for a successful fit.
%
% e.suturina@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=ilpcs.m>

function [mxyz,chi,Ilm,pred_pcs,s_mxyz,s_chi,s_Ilm]=ilpcs(nxyz,expt_pcs,ranks,mguess)

% Check consistency
grumble(nxyz,expt_pcs,ranks,mguess);

% Count the multipole variables
n_mvars=sum(2*nonzeros(ranks)+1);

% Set optimisation options
warning('off','optim:fminunc:SwitchingMethod');
options=optimset('Display','iter','MaxIter',inf,'UseParallel',true,...
                 'MaxFunEvals',inf,'FinDiffType','central');

% Vector of squared residuals 
vec_res_sq=@(x)(expt_pcs-lpcs(nxyz,x(1:3),ranks,...
                multipack(ranks,[0.5/sqrt(pi) x(9:end)]),x(4:8))).^2;
             
% Sum of squared residuals 
sum_res_sq=@(x)sum(vec_res_sq(x));

% Run the fitting to the multipole model
p=fminunc(sum_res_sq,[mguess 0.1 0.1 0.1 0.1 0.1 zeros(1,n_mvars)],options);

% Get metal coordinates
mxyz=p(1:3);

% Get susceptibility
chi=[p(4)     p(5)     p(6); 
     p(5)     p(7)     p(8); 
     p(6)     p(8) -p(4)-p(7)];
 
% Get multipoles
Ilm=multipack(ranks,[0.5/sqrt(pi) p(9:end)]);

% Get the predicted PCS values
pred_pcs=lpcs(nxyz,p(1:3),ranks,Ilm,p(4:8));

% If necessary, compute the statistics
if nargout>4

    % Compute Jacobian at the optimal point
    jac=jacobianest(vec_res_sq,p);

    % Get the Studentized residual
    sdr=sqrt(sum_res_sq(p)/(numel(expt_pcs)-n_mvars-8));

    % Get the standard deviations
    sp=sqrt(diag((sdr^2)*inv(jac'*jac)))'; %#ok<MINV>

    % Get metal SD
    s_mxyz=sp(1:3); 

    % Get susceptibility SD
    s_chi=[sp(4)          sp(5)          sp(6);
           sp(5)          sp(7)          sp(8);
           sp(6)          sp(8)  sqrt(sp(4)^2+sp(7)^2)];
       
    % Get multipole SD
    s_Ilm=multipack(ranks,[0 sp(9:end)]);

end

end

% Consistency enforcement
function grumble(nxyz,expt_pcs,L,mguess)
if (~isnumeric(nxyz))||(~isreal(nxyz))||(size(nxyz,2)~=3)
    error('nxyz must be an Nx3 array of atomic coordinates.');
end
if (~isnumeric(expt_pcs))||(size(expt_pcs,2)~=1)||(~isreal(expt_pcs))
    error('pcs_vals parameter should be a real column vector.');
end
if size(nxyz,1)~=size(expt_pcs,1)
    error('the number of rows in nxyz and pcs_vals arguments must be the same.');
end
if (~isnumeric(mguess))||(size(mguess,2)~=3)||(size(mguess,1)~=1)||(~isreal(mguess))
    error('mxyz parameter should be a real row vector with three elements.');
end
if (~isnumeric(L))||(~isreal(L))
    error('L must be a real vector.');
end
end

% As for butter versus margarine, I trust cows
% more than chemists.
%
% Joan Gussow

