% Fourier spectral differentiation as three polyadic factors. Syntax:
%
%                    D=fourdif_poly(npoints,order,extent)
%
% Parameters:
%
%    npoints - number of equally spaced periodic grid points
%    order   - positive integer derivative order
%    extent  - positive period length in the coordinate units
%
% Outputs:
%
%    D       - implicit polyadic derivative on column vectors or blocks
%
% The action is inverse FFT times a numeric CPU diagonal multiplier
% times FFT. Adjoint normalisation follows MATLAB conventions. Even
% grids use the zero Nyquist convention for odd derivatives, matching
% fourdif; even derivatives retain that mode. Upload D once with
% gpuArray when GPU execution is required. Do not inflate D.
%
% ilya.kuprov@weizmann.ac.il

function D=fourdif_poly(npoints,order,extent)

% Check consistency
grumble(npoints,order,extent);

% Build the spectral multiplier once on CPU
mult=fftdiff(order,npoints,extent/npoints).';
if (mod(npoints,2)==0)&&(mod(order,2)==1)
    mult(npoints/2+1)=0;
end
mult=spdiags(mult,0,npoints,npoints);

% Define transforms and their normalised adjoints
fft_core=struct('action',@(x)fft(x,[],1),...
                'adjoint',@(x)npoints*ifft(x,[],1),'dims',[npoints npoints]);
ifft_core=struct('action',@(x)ifft(x,[],1),...
                 'adjoint',@(x)fft(x,[],1)/npoints,'dims',[npoints npoints]);

% Compose inverse FFT, multiplication, and FFT factors
D=polyadic({{ifft_core}})*polyadic({{mult}})*polyadic({{fft_core}});

end

% Consistency enforcement
function grumble(npoints,order,extent)
if (~isnumeric(npoints))||(~isreal(npoints))||(~isscalar(npoints))||...
   (~isfinite(npoints))||(npoints<1)||(mod(npoints,1)~=0)
    error('npoints must be a finite positive integer.');
end
if (~isnumeric(order))||(~isreal(order))||(~isscalar(order))||...
   (~isfinite(order))||(order<1)||(mod(order,1)~=0)
    error('order must be a finite positive integer.');
end
if (~isnumeric(extent))||(~isreal(extent))||(~isscalar(extent))||...
   (~isfinite(extent))||(extent<=0)
    error('extent must be a finite positive real scalar.');
end
end


