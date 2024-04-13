% Basic 1D imaging sequence with a hard pulse. Syntax:
%
%       fid=basic_1d_hard(spin_system,parameters,H,R,K,G,F)
%
% This sequence must be called from the imaging() context, which
% would provide H,R,K,G, and F. Parameters:
%
%    parameters.ro_grad_amp  - readout gradient amplitude, T/m
%
%    parameters.sweep        - detection sweep width, Hz
%
%    parameters.npoints      - number of points in the fid
%
%    parameters.offset       - transmitter and receiver offset, Hz
%
% Outputs:
%
%    fid - free induction decay that should be Fourier transformed
%          to obtain the image
%
% i.kuprov@soton.ac.uk
% a.j.allami@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=basic_1d_hard.m>

function fid=basic_1d_hard(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,G,F);

% Assemble the Liouvillian
L=H+F+1i*R+1i*K;

% Make pulse operators
Hp=operator(spin_system,'L+',parameters.spins{1});
Hy=kron(speye(parameters.npts),(Hp-Hp')/2i);
Hx=kron(speye(parameters.npts),(Hp+Hp')/2);

% Hard 90-degree pulse
parameters.rho0=step(spin_system,Hy,parameters.rho0,pi/2);

% Echo time
parameters.rho0=evolution(spin_system,L,[],parameters.rho0,...
                          0.5/parameters.sweep,parameters.npoints,'final');
                      
% Hard 180-degree pulse
parameters.rho0=step(spin_system,Hx,parameters.rho0,pi);

% Pre-phasing gradient
parameters.rho0=evolution(spin_system,L-parameters.ro_grad_amp*G{1},[],...
                          parameters.rho0,0.5/parameters.sweep,parameters.npoints,'final');

% Acquisition under a gradient
fid=acquire(spin_system,parameters,L+parameters.ro_grad_amp*G{1},sparse(0),sparse(0));
         
end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K,G,F) %#ok<INUSL>
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('this function is only available in sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~isnumeric(F))||(~ismatrix(H))||(~ismatrix(R))||...
   (~ismatrix(K))||(~ismatrix(F))
    error('H,R,K,F arguments must be matrices.');
end
if (~all(size(H)==size(R)))||...
   (~all(size(R)==size(K)))||...
   (~all(size(K)==size(F)))
    error('H,R,K,F matrices must have the same dimension.');
end
if ~iscell(G)
    error('the G argument must be a cell array.');
end
end

% The term "Cobra Effect" stems from an anecdote set at the time
% of British rule of colonial India. The British government was
% concerned about the number of venomous cobra snakes in Delhi.
% The government therefore offered a bounty for every dead cobra.
% Initially this was a successful strategy as large numbers of
% snakes were killed for the reward. Eventually, however, enter-
% prising people began to breed cobras for the income. When the
% government became aware of this, the reward program was scrap-
% ped, causing the cobra breeders to set the now-worthless snakes
% free. As a result, the wild cobra population further increased.
%
% Wikipedia

