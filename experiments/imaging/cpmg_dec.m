% The effect of Carr-Purcell-Meiboom-Gill (CPMG) pulse sequence
% on the MRI phantom. The function runs the CPMG and then pro-
% jects out the user-specified spin state, returning the corres-
% ponding image. Syntax: 
% 
%          mri=cpmg_dec(spin_system,parameters,H,R,K,G,F)
%
% This sequence must be called from the imaging() context, which
% would provide H,R,K,G, and F. Parameters:
%
%    parameters.dec_time     - total duration of the sequence
%
%    parameters.npulses      - number of pulses in the sequence,
%                              excluding the first pi/2 pulse
%
%    parameters.spins        - nuclei on which the sequence
%                              is to act, e.g. {'1H'}
%
% Outputs:
%
%    mri  - amplitude of the detection state at each point of the
%           sample
%
% Note: the spin state to be observed should be specified in
%       parameters.coil_st, the coil phantom is ignored.
%
% a.j.allami@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=cpmg_dec.m>

function mri=cpmg_dec(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,G,F);

% Assemble the background
B=H+F+1i*R+1i*K;

% Make pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i);
Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2);

% Get CPMG delay
cpmg_delay=parameters.dec_time/parameters.npulses;

% Apply the first pulse
rho=step(spin_system,Ly,parameters.rho0,pi/2);

% Apply the first delay
rho=step(spin_system,B,rho,cpmg_delay/2);

% Run CPMG echo train
for n=1:(parameters.npulses-1)
    
    % Apply the pulse
    rho=step(spin_system,Lx,rho,pi);
    
    % Apply the delay
    rho=step(spin_system,B,rho,cpmg_delay);
    
    % Inform the user
    report(spin_system,['CPMG pi pulse ' num2str(n) ' out of ' ...
                         num2str(parameters.npulses) ' done.']);
    
end

% Apply the last pulse
rho=step(spin_system,Lx,rho,pi);

% Apply the last delay
rho=step(spin_system,B,rho,cpmg_delay/2);

% Project out the spin state and get its phantom
mri=fpl2phan(rho,parameters.coil_st{1},parameters.npts);

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K,G,F)
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
if ~iscell(G), error('the G argument must be a cell array.'); end
if ~isfield(parameters,'dec_time')
    error('sequence duration must be specfied in parameters.dec_time field.');
end
if (~isreal(parameters.dec_time))||(~isscalar(parameters.dec_time))||...
   (parameters.dec_time<=0)
    error('parameters.dec_time must be a positive real number.');
end
if ~isfield(parameters,'npulses')
    error('number of pi pulses must be specfied in parameters.npulses field.');
end
if (~isreal(parameters.npulses))||(~isscalar(parameters.npulses))||...
   (parameters.npulses<1)||(mod(parameters.npulses,1)~=0)
    error('parameters.npulses must be a positive real integer.');
end
end

% Both authors contributed equally to this work. Order of 
% authorship was determined by reproductive fitness.
%
% Footnote in https://doi.org/10.1016/j.jtbi.2007.03.015

