% Performs free induction decay apodisation. Supports free induction decays
% of any dimension. To satisfy Fourier transform symmetry requirements, the
% first elements of the FID in each dimension are divided by 2, except for
% singleton dimensions and those the user designates inactive. Syntax:
%
%                 fid=apodisation(spin_system,fid,winfuns)
%
% Arguments:
%
%     fid      - the free induction decay. The function expects a column
%                vector in the case of 1D FID, a 2D matrix with the time
%                origin located at the (1,1) corner point in the case of
%                a 2D FID, a 3D matrix with the time origin located at 
%                the (1,1,1) corner point in the case of a 3D FID, etc.
%
%     winfuns  - a cell array of window function specifications for each
%                dimension of the FID in the format {{spec},{spec},...},
%                omitting singleton dimensions. The following specifica-
%                tions are supported:
%
%         {}           - do nothing in this dimension
%
%         {'none'}     - no window function, but divide the first
%                        point by 2 to satisfy the Fourier trans-
%                        form symmetry requirement
%
%         {'crisp'}    - multiplied by cos(x)^8 half-bell. First 
%                        point has x=0, last point has x=pi/2.
%
%         {'exp',k}    - multiplied by exp(-k*x). First point has 
%                        x=0, last point has x=1.
%
%         {'gauss',k}  - multiplied by exp(-k*(x.^2)). First point 
%                        has x=0, last point has x=1.
%
%         {'cos'}      - multiplied by cos(x) half-bell. First po-
%                        int has x=0, last point has x=pi/2.
%
%         {'sin'}      - multiplied by sin(x) full bell. First po-
%                        int has x=0, last point has x=pi.
%
%         {'sqcos'}    - multiplied by cos(x).^2 half-bell. First 
%                        point has x=0, last point has x=pi/2.
%
%         {'sqsin'}    - multiplied by sin(x).^2 full bell. First 
%                        point has x=0, last point has x=pi.
%
%         {'kaiser',k} - multiplied by a Kaiser function with the 
%                        side lobe attenuation factor k. The peak 
%                        of the Kaiser function is in the middle 
%                        of the FID.
%
% Outputs:
%
%     fid   - apodised free induction decay
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=apodisation.m>

function fid=apodisation(spin_system,fid,winfuns)

% Check consistency
grumble(fid,winfuns);

% Find non-singleton dimensions
rel_dims=true(1,ndims(fid));
rel_dims(size(fid)<2)=false();
rel_dims=find(rel_dims);

% Exclude inactive dimensions
inact_dims=find(cellfun(@isempty,winfuns));
rel_dims=setdiff(rel_dims,inact_dims);

% Factors of 2
for dim=rel_dims

    % Index first points along each dimension
    idx=repmat({':'},1,ndims(fid)); idx{dim}=1;

    % FFT symmetry requirement
    fid(idx{:})=fid(idx{:})/2;

    % Report to the user
    report(spin_system,['FID dimension ' num2str(dim) ...
                        ', all first points divided by 2']);

end

% Window functions
for n=1:numel(rel_dims)

    % Get dimension and its point count
    dim=rel_dims(n); npts=size(fid,dim);

    % Build window function
    switch winfuns{dim}{1}

        case 'none'

            wf=ones(npts,1);

        case 'crisp'

            x=linspace(0,pi/2,npts);
            wf=cos(x(:)).^8;

        case 'exp'

            x=linspace(0,1,npts);
            k=winfuns{dim}{2};
            wf=exp(-k*x(:));

        case 'gauss'

            x=linspace(0,1,npts);
            k=winfuns{dim}{2};
            wf=exp(-k*(x(:).^2));

        case 'cos'

            x=linspace(0,pi/2,npts);
            wf=cos(x(:));

        case 'sin'

            x=linspace(0,pi,npts);
            wf=sin(x(:));

        case 'sqcos'

            x=linspace(0,pi/2,npts);
            wf=cos(x(:)).^2;

        case 'sqsin'

            x=linspace(0,pi,npts);
            wf=sin(x(:)).^2;

        case 'kaiser'

            k=winfuns{dim}{2};
            wf=kaiser(npts,k);

        otherwise

            % Complain and bomb out
            error(['window function type ' winfuns{dim}{1} ...
                ' is not implemented.']);

    end

    % Apply in the active dimension
    idx=ones(1,ndims(fid)); idx(dim)=npts;
    fid=fid.*reshape(wf,idx);

    % Report to the user
    report(spin_system,['FID dimension ' num2str(dim) ', '...
                        winfuns{dim}{1} ' window function applied.']);

end

end

% Consistency enforcement
function grumble(fid,winfuns)
if ~isnumeric(fid)
    error('fid must be a numeric array.');
end
if ~iscell(winfuns)
    error('winfuns must be a cell array.');
end
for n=1:numel(winfuns)
    if ~iscell(winfuns{n})
        error('elements of winfuns must be cell arrays.');
    end
end
end

% I have had my results for a long time, but I do not yet
% know how I am to arrive at them.
% 
% Carl Friedrich Gauss

