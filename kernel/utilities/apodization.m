% Performs free induction decay apodization. Supports 1D, 2D, and 3D free
% induction decays. Syntax:
%
%                     fid=apodization(fid,window_type,k)
%
% Arguments:
%
%     fid         - the free induction decay. The function expects a column
%                   vector in the case of 1D FID, a 2D matrix with the time
%                   origin located at the (1,1) corner point in the case of
%                   a 2D FID, and a 3D matrix with the time origin located
%                   at the (1,1,1) corner point in the case of a 3D FID.
%
%     window_type - type of the window function. The following window func-
%                   tion types are supported:
%
% One-dimensional FIDs and their arrays (time is the first dimension, array
% is over the second dimension). The first point of each FID in the array is
% divided by 2, and then each FID is:
%
%           'none-1d' - left unchanged
%
%          'crisp-1d' - multiplied by by cos(x)^8 half-bell. First 
%                       point has x=0, last point has x=pi/2.
%
%            'exp-1d' - multiplied by exp(-k*x). First point has x=0,
%                       last point has x=1.
%
%       'gaussian-1d' - multiplied by exp(-k*(x.^2)). First point has 
%                       x=0, last point has x=1.
%
%        'cosbell-1d' - multiplied by cos(x) half-bell. First point
%                       has x=0, last point has x=pi/2.
%
%        'sinbell-1d' - divides the first point by 2 and multiplies
%                       the FID by sin(x) full bell. First point
%                       has x=0, last point has x=pi.
%
%      'sqcosbell-1d' - multiplied by by cos(x).^2 half-bell. First 
%                       point has x=0, last point has x=pi/2.
%
%      'sqsinbell-1d' - multiplied by sin(x).^2 full bell. First point
%                       has x=0, last point has x=pi.
%
%         'kaiser-1d' - multiplied by a Kaiser function with the side
%                       lobe attenuation factor k. The peak of the
%                       Kaiser function is in the middle of the fid.
%
% Two-dimensional FIDs. The first column is divided by 2, then the first row
% is divided by 2, and then each FID is:
%
%           'none-2d' - left unchanged
%
%          'crisp-2d' - multiplied by cos(x)^4 half-bell in both 
%                       dimensions. First point in each dimension
%                       has x=0, last point has x=pi/2.
%
%            'exp-2d' - multiplied by exp(-k*x) in both dimensions.
%                       First point in each dimension has x=0, last
%                       point has x=1.
%
%       'gaussian-2d' - multiplied by exp(-k*(x.^2)) in both dimen-
%                       sions. First point in each dimension has
%                       x=0, last point has x=1.
%
%        'cosbell-2d' - multiplied by cos(x) half-bell in both dimen-
%                       sions. First point in each dimension has x=0,
%                       last point has x=pi/2.
%
%      'sqcosbell-2d' - multiplied by cos(x).^2 half-bell in both 
%                       dimensions. First point in each dimension 
%                       has x=0, last point has x=pi/2.
%
%        'sinbell-2d' - multiplied by sin(x) full bell in both dimen-
%                       sions. First point in each dimension has x=0,
%                       last point has x=pi.
%
%      'sqsinbell-2d' - multiplied by sin(x).^2 full bell in both
%                       dimensions. First point in each dimension
%                       has x=0, last point has x=pi.
%
% Three-dimensional FIDs. The first column is divided by 2, then the first row
% is divided by 2, then the first stack is divided by 2, and then each FID is:
%
%        'cosbell-3d' - multiplied by cos(x) half-bell in all dimen-
%                       sions. First point in each dimension has x=0,
%                       last point has x=pi/2.
%
%      'sqcosbell-3d' - multiplied by cos(x).^2 half-bell in all di-
%                       mensions. First point in each dimension has
%                       x=0, last point has x=pi/2.
%
%     k - decay rate parameters for those window functions that require
%         such parameters. The units are the reciprocal duration of the
%         fid: effectively the last point of the fid is taken to be t=1.
%         If a window function does not require a parameter, the value
%         of k is ignored. Default value of k is zero.
%
% Outputs:
%
%     fid   - apodised free induction decay
%
% Note: if an array of 1D FIDs is supplied, the FIDs must be in columns.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=apodization.m>

function fid=apodization(fid,window_type,k)

% Set the defaults
if ~exist('k','var'), k=0; end

% Check consistency
grumble(fid,window_type,k);

% Apply the window function
switch window_type
    
    case 'none-1d'
        fid(1,:)=fid(1,:)/2;
    
    case 'crisp-1d'
        fid(1,:)=fid(1,:)/2;
        x=linspace(0,pi/2,size(fid,1))';
        fid=fid.*cos(x).^8;
    
    case 'exp-1d'
        fid(1,:)=fid(1,:)/2;
        x=linspace(0,1,size(fid,1))';
        fid=fid.*exp(-k*x);
    
    case 'gaussian-1d'
        fid(1,:)=fid(1,:)/2;
        x=linspace(0,1,size(fid,1))';
        fid=fid.*exp(-k*(x.^2));
        
    case 'cosbell-1d'
        fid(1,:)=fid(1,:)/2;
        x=linspace(0,pi/2,size(fid,1))';
        fid=fid.*cos(x);
        
    case 'sinbell-1d'
        fid(1,:)=fid(1,:)/2;
        x=linspace(0,pi,size(fid,1))';
        fid=fid.*sin(x);
        
    case 'sqcosbell-1d'
        fid(1,:)=fid(1,:)/2;
        x=linspace(0,pi/2,size(fid,1))';
        fid=fid.*cos(x).^2;
        
    case 'sqsinbell-1d'
        fid(1,:)=fid(1,:)/2;
        x=linspace(0,pi,size(fid,1))';
        fid=fid.*sin(x).^2;
        
    case 'kaiser-1d'
        fid(1,:)=fid(1,:)/2;
        fid=fid.*kaiser(size(fid,1),k);
        
    case 'none-2d'
        fid(1,:)=fid(1,:)/2;
        fid(:,1)=fid(:,1)/2;
        
    case 'crisp-2d'
        fid(1,:)=fid(1,:)/2;
        fid(:,1)=fid(:,1)/2;
        x=linspace(0,pi/2,size(fid,1))';
        y=linspace(0,pi/2,size(fid,2))';
        decay_col=cos(x).^4;
        decay_row=cos(y).^4;
        fid=fid.*kron(decay_col,decay_row');
    
    case 'exp-2d'
        fid(1,:)=fid(1,:)/2;
        fid(:,1)=fid(:,1)/2;
        x=linspace(0,1,size(fid,1))';
        y=linspace(0,1,size(fid,2))';
        decay_col=exp(-k*x);
        decay_row=exp(-k*y);
        fid=fid.*kron(decay_col,decay_row');
    
    case 'gaussian-2d'
        fid(1,:)=fid(1,:)/2;
        fid(:,1)=fid(:,1)/2;
        x=linspace(0,1,size(fid,1))';
        y=linspace(0,1,size(fid,2))';
        decay_col=exp(-k*(x.^2));
        decay_row=exp(-k*(y.^2));
        fid=fid.*kron(decay_col,decay_row');
        
    case 'cosbell-2d'
        fid(1,:)=fid(1,:)/2;
        fid(:,1)=fid(:,1)/2;
        x=linspace(0,pi/2,size(fid,1))';
        y=linspace(0,pi/2,size(fid,2))';
        decay_col=cos(x);
        decay_row=cos(y);
        fid=fid.*kron(decay_col,decay_row');
        
    case 'sqcosbell-2d'
        fid(1,:)=fid(1,:)/2;
        fid(:,1)=fid(:,1)/2;
        x=linspace(0,pi/2,size(fid,1))';
        y=linspace(0,pi/2,size(fid,2))';
        decay_col=cos(x).^2;
        decay_row=cos(y).^2;
        fid=fid.*kron(decay_col,decay_row');
    
    case 'sinbell-2d'
        fid(1,:)=fid(1,:)/2;
        fid(:,1)=fid(:,1)/2;
        x=linspace(0,pi,size(fid,1))';
        y=linspace(0,pi,size(fid,2))';
        decay_col=sin(x);
        decay_row=sin(y);
        fid=fid.*kron(decay_col,decay_row');
        
    case 'sqsinbell-2d'
        fid(1,:)=fid(1,:)/2;
        fid(:,1)=fid(:,1)/2;
        x=linspace(0,pi,size(fid,1))';
        y=linspace(0,pi,size(fid,2))';
        decay_col=sin(x).^2;
        decay_row=sin(y).^2;
        fid=fid.*kron(decay_col,decay_row');
        
    case 'cosbell-3d'
        fid(1,:,:)=fid(1,:,:)/2;
        fid(:,1,:)=fid(:,1,:)/2;
        fid(:,:,1)=fid(:,:,1)/2;
        [f1,f2,f3]=ndgrid(linspace(0,pi/2,size(fid,1)),...
                          linspace(0,pi/2,size(fid,2)),...
                          linspace(0,pi/2,size(fid,3)));
        fid=fid.*cos(f1).*cos(f2).*cos(f3);
        
    case 'sqcosbell-3d'
        fid(1,:,:)=fid(1,:,:)/2;
        fid(:,1,:)=fid(:,1,:)/2;
        fid(:,:,1)=fid(:,:,1)/2;
        [f1,f2,f3]=ndgrid(linspace(0,pi/2,size(fid,1)),...
                          linspace(0,pi/2,size(fid,2)),...
                          linspace(0,pi/2,size(fid,3)));
        fid=fid.*(cos(f1).^2).*(cos(f2).^2).*(cos(f3).^2);
        
    otherwise
        
        % Complain and bomb out
        error(['function ' window_type ' not implemented.']);

end

end

% Consistency enforcement
function grumble(fid,window_type,decay_rate)
if (~isnumeric(fid))||(~isnumeric(decay_rate))
    error('fid and decay rate must be numeric.');
end
if (numel(decay_rate)~=1)||(~isreal(decay_rate))||(decay_rate<0)
    error('decay rate must be a positive real number.');
end
if ~ischar(window_type)
    error('window_type must be a character string.');
end
if strcmp(window_type((end-1):end),'2d')&&...
   ((numel(size(fid))~=2)||any(size(fid)<2))
    error('2D window functions require FID to be a 2D matrix.');
end
if strcmp(window_type((end-1):end),'3d')&&...
   ((numel(size(fid))~=3)||any(size(fid)<2))
    error('3D window functions require FID to be a 3D matrix.');
end
end

% I have had my results for a long time, but I do not yet
% know how I am to arrive at them.
% 
% Carl Friedrich Gauss

