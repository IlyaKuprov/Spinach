% Computes rotational correlation functions using a Monte-Carlo method and
% compares them to the analytical results returned by Spinach kernel for
% the following correlation function:
%
%                    G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'>
%
% The sigma parameters refer to the rates of rotation and the four indices
% to the Wigner functions being correlated.
%
% High-rank isotropic rotational diffusion tested here.
%
% ilya.kuprov@weizmann.ac.il

function correlation_function_4()

% Set testing parameters
sigma_iso=0.2; L=4; k=-1; m=2; p=-1; q=2;

% Convert indices from [-L,L] to [1,2*L+1]
k=L+1+k; m=L+1+m; p=L+1+p; q=L+1+q;

%% Numerical Monte-Carlo calculation

% Number of points and lags
npoints=1e5; nlags=100;

% Generate angle track
angles=randn(3,npoints);

% Preallocate and start DCM trajectory
DCMT=zeros([3 3 npoints]); DCM=eye(3);

% Loop over Monte-Carlo steps
for n=1:npoints

    % Store the step
    DCMT(:,:,n)=DCM;
    
    % Generate a random rotation
    R_gen=sigma_iso*[ 0  1  0; -1  0  0;  0  0  0]*angles(1,n)+...
          sigma_iso*[ 0  0  1;  0  0  0; -1  0  0]*angles(2,n)+...
          sigma_iso*[ 0  0  0;  0  0  1;  0 -1  0]*angles(3,n);
    DCM=DCM*expm(R_gen);
      
end

% Convert DCMs into Wigner functions
W=zeros([2*L+1, 2*L+1, npoints],'like',1i);
parfor n=1:size(DCMT,3)

    % Pull a DCM out
    DCM=squeeze(DCMT(:,:,n));

    % Convert into Euler angles
    [alp,bet,gam]=dcm2euler(DCM);

    % Convert into Wigner functions
    W(:,:,n)=wigner(L,alp,bet,gam);

end

% Get Monte-Carlo correlation function
[cf_mc,lags]=xcorr(squeeze(W(k,m,:)),...
                   squeeze(W(p,q,:)),nlags,'normalized');
cf_mc=(1/(2*L+1))*ifftshift(cf_mc); lags=ifftshift(lags);

%% Analytical Spinach calculation

% Create a dummy spin system
sys.magnet=0; sys.isotopes={'G'};
inter.relaxation={'redfield'};
inter.rlx_keep='labframe';
inter.equilibrium='zero';
inter.tau_c={1/(3*sigma_iso^2)};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get the analytical correlation function
[weights,rates]=corrfun(spin_system,L,k,m,p,q);
cf_an=zeros(1,nlags);
for k=1:numel(weights{1})
    cf_an=cf_an+weights{1}(k)*exp(rates{1}(k)*(0:(nlags-1)));
end

% Plotting
figure(); 
plot(lags(1:nlags),real(cf_mc(1:nlags)),'ro'); hold on;
plot(lags(1:nlags),cf_an,'b-'); xlim('tight'); kgrid; 
kylabel('correlation function'); kxlabel('lag, points');
klegend({'Monte-Carlo','Spinach'},'Location','NorthEast');

end

