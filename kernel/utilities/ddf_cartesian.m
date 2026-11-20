% Apply an isolated or periodic Cartesian dipolar field operator. Syntax:
%
%                  field=ddf_cartesian(magn,npts,symbol)
%
% Parameters:
%
%    magn   - prod(npts)-by-3 Cartesian magnetisation density in A/m
%    npts   - 1-by-3 Cartesian voxel counts, first index fastest
%    symbol - multiplier returned by ddf_bounded() or ddf_periodic()
%
% Outputs:
%
%    field  - prod(npts)-by-3 secular rotating-frame field in tesla

function field=ddf_cartesian(magn,npts,symbol)

% Check consistency
grumble(magn,npts,symbol);

% Transform one Cartesian component at a time to limit padded storage
fft_size=size(symbol,1:3);
field=zeros(size(magn));
for n=1:3
    source=reshape(magn(:,n),npts);
    density=fftn(source,fft_size);
    local=ifftn(symbol.*density,'symmetric');
    local=local(1:npts(1),1:npts(2),1:npts(3));
    field(:,n)=local(:);
end

% Convert to tesla and impose the transverse secular tensor components
field=1.25663706127e-6*field.*[-1/2 -1/2 1];

end

% Consistency enforcement
function grumble(magn,npts,symbol)
validateattributes(npts,{'double'},{'size',[1 3],'real','finite','integer','positive'});
validateattributes(magn,{'double'},{'size',[prod(npts) 3],'real','finite','nonsparse'});
if (ndims(symbol)>3)||(~isequal(size(symbol,1:3),npts)&&~isequal(size(symbol,1:3),2*npts-1))
    error('symbol dimensions must match the periodic or isolated grid.');
end
validateattributes(symbol,{'double'},{'real','finite','nonsparse','nonempty'});
end


