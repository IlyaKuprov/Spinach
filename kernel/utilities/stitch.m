% Stitching function for bidirectionally propagated 3D NMR pulse
% sequences. Propagate your initial condition forward to some mid-
% point, propagate your detection state backward to the same mid-
% point and use this function to obtain the 3D free induction de-
% cay (http://dx.doi.org/10.1016/j.jmr.2014.04.002) at the price
% of two 2D simulations. Syntax:
%
%   fid=stitch(spin_system,L,rho_stack,coil_stack,...
%                                  mtp_oper,mtp_time,t1,t2,t3)
%
% Parameters:
%
%             L - spin system Liouvillian
% 
%     rho_stack - state vector stack from the forward part of
%                 the simulation
%
%    coil_stack - coil vector stack from the backward part of
%                 the sumulation
%
%      mec_oper - cell array of operators in the midpoint event
%                 chain, e.g. {Lx,L,Sy}
%
%      mec_time - cell array of durations of each event at the
%                 midpoint of the t2 evolution period
%
%     t1.nsteps - number of time steps in t1
%
%     t2.nsteps - number of time steps in t2
%
%     t3.nsteps - number of time steps in t3
%
%          tdir - time direction for state and coil propagation,
%                 the default is '+-'
%
% Outputs:
%
%           fid - three-dimensional free induction decay
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=stitch.m>

function fid=stitch(spin_system,L,rho_stack,coil_stack,...
                    mec_oper,mec_time,t1,t2,t3,tdir)

% Set default time directions
if ~exist('tdir','var'), tdir='+-'; end

% Check consistency
grumble(L,rho_stack,coil_stack,mec_oper,mec_time,t1,t2,t3);

% Preallocate the fid
fid=zeros(t3.nsteps,t2.nsteps,t1.nsteps,'like',1i);

% Compute time step propagators
P=propagator(spin_system,L,t2.timestep/2); Pct=P';

% Compute midpoint event propagator
Pm=speye(size(L));
for n=1:numel(mec_oper)
    Pm=propagator(spin_system,mec_oper{n},mec_time{n})*Pm;
end
Pm=clean_up(spin_system,Pm,spin_system.tols.prop_chop);

% Decide the device and style
if ismember('gpu',spin_system.sys.enable)
    
    % Upload the problem to GPU
    P=gpuArray(P); 
    Pm=gpuArray(Pm); 
    Pct=gpuArray(Pct); 
    rho_stack=gpuArray(rho_stack); 
    coil_stack=gpuArray(coil_stack);
    
    % Inform the user
    report(spin_system,'stitching will be done on GPU.');
    
else
    
    % Inform the user
    report(spin_system,'stitching will be done on CPU.');
    
end

% Stitch the trajectories
for k=1:t2.nsteps
    
    % Inform the user
    report(spin_system,['stitching forward and backward trajectories, step '...
                         num2str(k) '/' num2str(t2.nsteps) '...']);
       
    % Scalar product with the midpoint propagator
    fid(:,k,:)=gather(coil_stack'*(Pm*rho_stack));
    
    % Time directions
    switch tdir
        
        case '++'
            
            % Forward time evolution for rho
            rho_stack=P*rho_stack;
 
            % Forward time evolution for coil
            coil_stack=P*coil_stack;
            
        case '+-'
            
            % Forward time evolution for rho
            rho_stack=P*rho_stack;
            
            % Backward time evolution for coil
            coil_stack=Pct*coil_stack;
            
        case '-+'
            
            % Backward time evolution for rho
            rho_stack=Pct*rho_stack;
            
            % Forward time evolution for coil
            coil_stack=P*coil_stack;
            
        case '--'
            
            % Backward time evolution for rho
            rho_stack=Pct*rho_stack;
            
            % Backward time evolution for coil
            coil_stack=Pct*coil_stack;
            
        otherwise
            
            % Complain and bomb out
            error('invalid time direction specification in tdir');
            
    end
    
end

end

% Consistency enforcement
function grumble(L,rho_stack,coil_stack,mec_oper,mec_time,t1,t2,t3)
if (~isnumeric(L))||(~ismatrix(L))||(size(L,1)~=size(L,2))
    error('L must be a square matrix.');
end
if (~isfield(t1,'nsteps'))||(~isnumeric(t1.nsteps))||...
   (~isreal(t1.nsteps))||(mod(t1.nsteps,1)~=0)||(t1.nsteps<1)
    error('t1.nsteps must exist and must be a positive real integer.');
end
if (~isfield(t2,'nsteps'))||(~isnumeric(t2.nsteps))||...
   (~isreal(t2.nsteps))||(mod(t2.nsteps,1)~=0)||(t2.nsteps<1)
    error('t2.nsteps must exist and must be a positive real integer.');
end
if (~isfield(t3,'nsteps'))||(~isnumeric(t3.nsteps))||...
   (~isreal(t3.nsteps))||(mod(t3.nsteps,1)~=0)||(t3.nsteps<1)
    error('t3.nsteps must exist and must be a positive real integer.');
end
if (~isnumeric(rho_stack))||(size(rho_stack,1)~=size(L,1))||...
   (size(rho_stack,2)~=t1.nsteps)
    error('rho_stack must be a numerical array with dimensions size(L,1) x t1.nsteps');
end
if (~isnumeric(coil_stack))||(size(coil_stack,1)~=size(L,1))||...
   (size(coil_stack,2)~=t3.nsteps)
    error('coil_stack must be a numerical array with dimensions size(L,1) x t3.nsteps');
end
if ~iscell(mec_oper)
    error('mtp_oper must be a cell array of matrices.');
end
for n=1:numel(mec_oper)
    if (~isnumeric(mec_oper{n}))||(~ismatrix(mec_oper{n}))||...
       (~all(size(mec_oper{n})==size(L)))
        error('mtp_oper contain square matrices of dimensions same as L.');
    end
end
if ~iscell(mec_time)
    error('mec_time must be a cell array of scalars.');
end
for n=1:numel(mec_time)
    if (~isnumeric(mec_time{n}))||(~isreal(mec_time{n}))||...
       (~isscalar(mec_time{n}))||(mec_time{n}<=0)
        error('elements of mec_time must be a positive real scalars.');
    end
end
end

% A wolf hates both men and dogs, but dogs he hates more.
%
% Sergey Dovlatov

