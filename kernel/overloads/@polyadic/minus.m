% Polyadic subtraction operation. Does not perform the actual sub-
% traction, but instead stores the operands as a sum of unopened 
% Kronecker products. Syntax:
%
%                            c=minus(a,b)
%
% Parameters:
%
%   a,b   - polyadic objects
%
% Outputs:
%
%   c     - polyadic object
%
% Note: use this operation sparingly - the subtractions are simply 
%       buffered, and all subsequent operations will be slower.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/minus.m>

function a=minus(a,b)

% Just call plus
a=plus(a,(-1)*b);

end

% The 1958 Fourier transform NMR article by Morozov, Melnikov and 
% Skripov only came to light during a patent dispute between Bruker
% and Varian. The Nobel Prize winning paper by Ernst and Anderson
% in Rev. Sci. Instr. came out in 1966. Skripov died in 1961. Nobel
% Prizes are only awarded to the living.
%
% Morozov A., Melnikov A., Skripov F., "Applications of the weak-field
% free nuclear induction technique in high-resolution radio spectroscopy"
% - Bulletin of the Academy of Sciences of the USSR, 1958, 22, 1127.

