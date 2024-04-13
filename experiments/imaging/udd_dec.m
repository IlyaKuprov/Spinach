% The effect of Uhrid Dynamic Decoupling (UDD) pulse sequence
% on the MRI phantom. The function runs the UDD and then pro-
% jects out the user-specified spin state, returning the cor-
% responding image. Syntax: 
% 
%          mri=udd_dec(spin_system,parameters,H,R,K,G,F)
%
% This sequence must be called from the imaging() context, which
% would provide H,R,K,G, and F. Parameters:
%
%    parameters.dec_time     - total duration of the sequence
%
%    parameters.npulses      - number of pulses in the sequence,
%                              including the first pi/2 pulse
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
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=udd_dec.m>

function mri=udd_dec(spin_system,parameters,H,R,K,G,F)

% Check consistency
grumble(spin_system,parameters,H,R,K,G,F);

% Assemble the background
B=H+F+1i*R+1i*K;

% Make pulse operators
Lp=operator(spin_system,'L+',parameters.spins{1});
Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i);
Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2);

% Apply 90-degree pulse
rho=step(spin_system,Ly,parameters.rho0,pi/2);

% Get UDD delays
udd_delays=uhrig_times(parameters.dec_time,parameters.npulses);

% Run UDD echo train
for n=1:(numel(udd_delays)-1)
    
    % Apply the delay
    rho=step(spin_system,B,rho,udd_delays(n));
    
    % Apply the pulse
    rho=step(spin_system,Lx,rho,pi);
    
    % Inform the user
    report(spin_system,['UDD pi pulse ' num2str(n) ' out of ' ...
                        num2str(numel(udd_delays)-1) ' done.']); 

end

% Apply the last delay
rho=step(spin_system,B,rho,udd_delays(end));

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

% Of course, the ideas of different members of the University community
% will often and quite naturally conflict. But it is not the proper role
% of the University to attempt to shield individuals from ideas and opi-
% nions they find unwelcome, disagreeable, or even deeply offensive. Al-
% though the University greatly values civility, and although all members
% of the University community share in the responsibility for maintaining
% a climate of mutual respect, concerns about civility and mutual respect
% can never be used as a justification for closing off discussion of ideas,
% however offensive or disagreeable those ideas may be to some members of
% our community.
%
% University of Chicago, on freedom of expression

