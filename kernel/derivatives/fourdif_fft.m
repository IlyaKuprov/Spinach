% Matrix-free first Fourier derivative on a periodic angular grid.
% Syntax:
%
%                         D=fourdif_fft(npoints)
%
% Parameters:
%
%    npoints - positive integer number of equispaced rotor grid points
%
% Outputs:
%
%    D       - matfree core representing fourdif(npoints,1)'s derivative;
%              a scalar zero for the one-point grid
%
% The action accepts full npoints-by-n blocks, real or complex, on CPU
% or GPU. Even grids use zero Nyquist derivative, matching fourdif.
% The angular coordinate spans 2*pi, so no grid-spacing factor is needed.
%
% ilya.kuprov@weizmann.ac.il

function D=fourdif_fft(npoints)

% Check consistency
grumble(npoints);

% A single rotor phase has zero derivative
if npoints==1, D=sparse(1,1); return; end

% Reuse the spectral kernel with the first-derivative Nyquist convention
kernel=fftdiff(1,npoints,2*pi/npoints).';
if mod(npoints,2)==0, kernel(npoints/2+1)=0; end

% The first derivative is real and skew-adjoint
forward=@(block)ifft(kernel.*fft(block,[],1),[],1);
adjoint=@(block)ifft(conj(kernel).*fft(block,[],1),[],1);
D=matfree([npoints npoints],forward,adjoint,true);

end

% Consistency enforcement
function grumble(npoints)
validateattributes(npoints,{'double'},{'scalar','integer','positive','finite'});
end


