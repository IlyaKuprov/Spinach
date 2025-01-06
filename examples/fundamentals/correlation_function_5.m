% Computes the following rotational correlation function
%
%                G(k,m,p,q)=<R(k,m)*R(p,q)>
%
% where R is the 3D Cartesian rotation matrix, using the 
% Monte-Carlo method.
%
% ilya.kuprov@weizmann.ac.il

function correlation_function_5()

% Set testing parameters
sigma_iso=0.2; k=2; m=3; p=2; q=3;

% Set number of points
npoints=1e5; nlags=300;

% Generate angle track
angles=randn(3,npoints);

% Preallocate rotation matrix array
R=zeros([3 3 npoints]); R(:,:,1)=eye(3);

% Loop over Monte-Carlo steps
for n=2:npoints
    
    % Generate a random rotation
    R_gen=[ 0  1  0; -1  0  0;  0  0  0]*angles(1,n)+...
          [ 0  0  1;  0  0  0; -1  0  0]*angles(2,n)+...
          [ 0  0  0;  0  0  1;  0 -1  0]*angles(3,n);
    R(:,:,n)=R(:,:,n-1)*expm(sigma_iso*R_gen);
       
end

% Get Monte-Carlo correlation function
[cf_mc,lags]=xcorr(squeeze(R(k,m,:)),...
                   squeeze(R(p,q,:)),nlags,'normalized');
cf_mc=(1/3)*ifftshift(cf_mc); lags=ifftshift(lags);

% Plotting
figure(); 
plot(lags(1:nlags),real(cf_mc(1:nlags)),'b-');
xlim('tight'); kgrid; kxlabel('lag, points');
kylabel('correlation function'); 
klegend({'Monte-Carlo'},'Location','NorthEast');

end

