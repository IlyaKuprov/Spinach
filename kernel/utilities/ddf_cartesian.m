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
if (~isa(npts,'double'))||(~isequal(size(npts),[1 3]))||(~isreal(npts))||...
   (any(~isfinite(npts(:))))||(any(npts(:)~=fix(npts(:))))||(any(npts(:)<=0))
    error('npts must be a 1-by-3 double row of positive integers.');
end
if (~isa(magn,'double'))||(~isequal(size(magn),[prod(npts) 3]))||(~isreal(magn))||...
   (any(~isfinite(magn(:))))||(issparse(magn))
    error('magn must be a full real double array with three finite components per voxel.');
end
if (ndims(symbol)>3)||(~isequal(size(symbol,1:3),npts)&&~isequal(size(symbol,1:3),2*npts-1))
    error('symbol dimensions must match the periodic or isolated grid.');
end
if (~isa(symbol,'double'))||(~isreal(symbol))||(any(~isfinite(symbol(:))))||(issparse(symbol))||...
   (isempty(symbol))
    error('symbol must be a real finite full non-empty double array.');
end
end


