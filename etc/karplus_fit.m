% Fits a Karplus curve to a Gaussian dihedral angle scan. Syntax:
%
%           [A,B,C,sA,sB,sC]=karplus_fit(dir_path,atoms)
%
% Parameters:
%
%    dir_path - path to the directory containing the 
%               Gaussian logs
%
%    atoms    - a cell array of 4-element vectors
%               specifying atmos making up the dihe-
%               dral angles of interest
%
% The directory specified in the first argument should contain 
% a series of Gaussian J-coupling calculation logs that differ 
% only in the value of the dihedral angle in question.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=Karplus_fit.m>

function [A,B,C,sA,sB,sC]=karplus_fit(dir_path,atoms)

% Check consistency
grumble(dir_path,atoms);

% Get all log files in the directory
logfiles=dir([dir_path '/*.log']);

% Preallocate the arrays
phi=zeros(numel(atoms),numel(logfiles));
J=zeros(numel(atoms),numel(logfiles));

% Parse the files
parfor n=1:numel(logfiles)
    props{n}=gparse([dir_path '/' logfiles(n).name]);
end

% Extract parameters
for n=1:numel(logfiles) %#ok<*PFBNS>
    try
        for k=1:numel(atoms)
            phi(k,n)=dihedral(props{n}.std_geom(atoms{k}(1),:),props{n}.std_geom(atoms{k}(2),:),...
                              props{n}.std_geom(atoms{k}(3),:),props{n}.std_geom(atoms{k}(4),:)); 
            J(k,n)=props{n}.j_couplings(atoms{k}(1),atoms{k}(4));
        end
    catch
        disp(['Gaussian log file ' logfiles(n).name ' could not be interpreted, skipped.']);
        for k=1:numel(atoms)
            phi(k,n)=NaN; J(k,n)=NaN;
        end
    end
end

% Eliminate NaNs
mask=isnan(phi)|isnan(J);
phi(mask)=[]; J(mask)=[]; 

% Stretch the arrays
phi=phi(:)'; J=J(:)';

% Rotate phi into [0,360]
phi=mod(phi,360);

% Compute basis cosines
cosine_2=cosd(phi).^2;
cosine_1=cosd(phi);
cosine_0=ones(size(phi));

% Run linear least squares
result=[cosine_2' cosine_1' cosine_0']\J';
A=result(1); B=result(2); C=result(3);

% Vector of squared residuals
vec_res_sq=@(x)(J'-x(1)*cosine_2'-x(2)*cosine_1'-x(3)*cosine_0').^2;

% Sum of squared residuals 
sum_res_sq=@(x)sum(vec_res_sq(x));

% Compute the Jacobian at the optimal point
jac=jacobianest(vec_res_sq,result);

% Get the Studentized residual
sdr=sqrt(sum_res_sq(result)/(numel(J)-3));

% Get the standard deviations
stdevs=sqrt(diag((sdr^2)*inv(jac'*jac))); %#ok<MINV>
sA=stdevs(1); sB=stdevs(2); sC=stdevs(3);

% Plot Karplus curve
figure(); plot(phi,J,'ro'); kgrid;
psi=linspace(0,360,128); hold on; 
plot(psi,A*cosd(psi).^2+B*cosd(psi)+C,'b-');
xlabel('Dihedral angle, degrees');
ylabel('J-coupling, Hz'); axis tight;

end

% Consistency enforcement
function grumble(dir_path,atoms)
if ~ischar(dir_path)
    error('dir_path must be a character string.');
end
if ~iscell(atoms)
    error('atoms must be a cell array of 4-element vectors.');
end
end

% If men learn this, it will implant forgetfulness in their souls; they will
% cease to exercise memory because they rely on that which is written [...]
% it is no true wisdom that you offer your disciples, but only its semblance,
% for by telling them of many things without teaching them you will make them
% seem to know much, while for the most part they know nothing, and as men 
% filled, not with wisdom but with the conceit of wisdom, they will be a bur-
% den to their fellows.
%
% Plato (circa 429-347 BCE), about reading and writing.

