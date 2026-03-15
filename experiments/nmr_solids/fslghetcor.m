% Heteronuclear correlation MAS NMR experiment with frequency-switched
% Lee-Goldburg homonuclear decoupling. Further details in:
% 
%    https://doi.org/10.1103/PhysRev.140.A1261
%    https://doi.org/10.1016/0009-2614(89)87166-0
%    https://doi.org/10.1006/jmre.1996.1089 (Figure 1)
%
% Syntax:
%
%              fid=fslghetcor(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.spins        - working spins, e.g. {'1H',13C'}
%
%     parameters.hi_pwr       - amplitude of high power pulses 
%                               on the high-gamma channel, Hz
%
%     parameters.cp_pwr       - amplitude of CP pulse on each 
%                               channel during the CP contact
%                               time, Hz
%
%     parameters.cp_dur       - CP contact time duration, s
%
%     parameters.offset       - transmitter offsets on the
%                               two channels, Hz
%
%     parameters.nblocks      - number of FSLG blocks per
%                               indirect-dimension point
%
%     parameters.spc_dim      - Fokker-Planck spatial dimension
%
%     parameters.rho0         - initial state
%
%     parameters.coil         - detection state
%
%     parameters.sweep        - sweep width, Hz for F1, F2
%
%     parameters.npoints      - number of points in F1, F2
%
%     H - Hamiltonian superoperator, received from context function
%
%     R - relaxation superoperator, received from context function
%
%     K - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid.sin, fid.cos        - sine and cosine components
%                              of the States quadrature
%
% guinevere.mathies@uni-konstanz.de
%
% <https://spindynamics.org/wiki/index.php?title=fslghetcor.m>

function fid=fslghetcor(spin_system,parameters,H,R,K)
 
% Consistency enforcement
grumble(spin_system,parameters,H,R,K);
 
% Build 1H and 13C control operators
Hx=operator(spin_system,'Lx',parameters.spins{1}); Hx=kron(speye(parameters.spc_dim),Hx);
Hy=operator(spin_system,'Ly',parameters.spins{1}); Hy=kron(speye(parameters.spc_dim),Hy);
Hz=operator(spin_system,'Lz',parameters.spins{1}); Hz=kron(speye(parameters.spc_dim),Hz);
Cx=operator(spin_system,'Lx',parameters.spins{2}); Cx=kron(speye(parameters.spc_dim),Cx);

% Liouvillian 
L=H+1i*R+1i*K;

% Proton pulse evolution generators
L_HPx=L+2*pi*parameters.hi_pwr*Hx; L_HPy=L+2*pi*parameters.hi_pwr*Hy;
L_HMx=L-2*pi*parameters.hi_pwr*Hx; L_HMy=L-2*pi*parameters.hi_pwr*Hy;

% FSLG pulse evolution generators
L_FSLG_p=L-2*pi*parameters.offset(1)*Hz ...
          +2*pi*(parameters.hi_pwr/sqrt(2))*Hz;
L_FSLG_m=L-2*pi*parameters.offset(1)*Hz ...
          -2*pi*(parameters.hi_pwr/sqrt(2))*Hz;
L1_cos=L_FSLG_m+2*pi*parameters.hi_pwr*Hy;
L2_cos=L_FSLG_p-2*pi*parameters.hi_pwr*Hy;
L1_sin=L_FSLG_p+2*pi*parameters.hi_pwr*Hx;
L2_sin=L_FSLG_m-2*pi*parameters.hi_pwr*Hx;

% CP period evolution generator
L_c=L-2*pi*parameters.cp_pwr(1)*Hy...
     +2*pi*parameters.cp_pwr(2)*Cx;

% Sequence timing shorthands
hi_pwr_90deg=1/(4*parameters.hi_pwr);
hi_pwr_magic=acos(1/sqrt(3))/(2*pi)/parameters.hi_pwr;
dwell_times=1./parameters.sweep;

% Flip by 90 degrees + magic angle
rho_cos=step(spin_system,L_HPx,parameters.rho0,hi_pwr_90deg+hi_pwr_magic);
rho_sin=step(spin_system,L_HPy,parameters.rho0,hi_pwr_90deg+hi_pwr_magic);

% Preallocate and start the F1 trajectories
traj_cos=zeros(numel(rho_cos),parameters.npoints(1),'like',1i);
traj_sin=zeros(numel(rho_sin),parameters.npoints(1),'like',1i);
traj_cos(:,1)=rho_cos; traj_sin(:,1)=rho_sin;

% Compute the F1 trajectory, cos part
if ismember('gpu',spin_system.sys.enable)
    L1_cos=gpuArray(L1_cos); L2_cos=gpuArray(L2_cos);
end
for k=2:parameters.npoints(1)
    for n=1:parameters.nblocks
        rho_cos=step(spin_system,L1_cos,rho_cos,sqrt(2/3)/parameters.hi_pwr);
        rho_cos=step(spin_system,L2_cos,rho_cos,sqrt(2/3)/parameters.hi_pwr);
    end
    traj_cos(:,k)=rho_cos;
end
clear('L1_cos','L2_cos');

% Compute the F1 trajectory, sin part
if ismember('gpu',spin_system.sys.enable)
    L1_sin=gpuArray(L1_sin); L2_sin=gpuArray(L2_sin);
end
for k=2:parameters.npoints(1)
    for n=1:parameters.nblocks
        rho_sin=step(spin_system,L1_sin,rho_sin,sqrt(2/3)/parameters.hi_pwr);
        rho_sin=step(spin_system,L2_sin,rho_sin,sqrt(2/3)/parameters.hi_pwr);
    end
    traj_sin(:,k)=rho_sin;
end
clear('L1_sin','L2_sin');

% Flip the 1H back from the magic angle to the x,y-plane
traj_cos=step(spin_system,L_HMx,traj_cos,hi_pwr_magic);
traj_sin=step(spin_system,L_HMy,traj_sin,hi_pwr_magic);

% Save work by stacking cos and sin parts
traj_cos_sin=[traj_cos traj_sin];

% Run the CP contact period
traj_cos_sin=step(spin_system,L_c,traj_cos_sin,parameters.cp_dur);
   
% Decouple 1H for the acquisition period
[L_dec,traj_cos_sin]=decouple(spin_system,L,traj_cos_sin,parameters.spins(1));
    
% Run the F2 evolution
fid_cos_sin=evolution(spin_system,L_dec,parameters.coil,traj_cos_sin,...
                      dwell_times(2),parameters.npoints(2)-1,'observable');

% Unstack cos and sin parts
fid.cos=fid_cos_sin(:,1:(end/2)); 
fid.sin=fid_cos_sin(:,(end/2+1):end); 

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'hi_pwr')
    error('high RF amplitude must be specified in parameters.hi_pwr variable.');
elseif (~isnumeric(parameters.hi_pwr))||(~isreal(parameters.hi_pwr))||...
       (numel(parameters.hi_pwr)~=1)||(parameters.hi_pwr<=0)
    error('parameters.hi_pwr must be a positive real scalar.');
end
if ~isfield(parameters,'cp_pwr')
    error('RF amplitude during CP must be specified in parameters.cp_pwr variable.');
elseif (~isnumeric(parameters.cp_pwr))||(~isreal(parameters.cp_pwr))||...
       (numel(parameters.cp_pwr)~=2)||any(parameters.cp_pwr<=0)
    error('parameters.cp_pwr must be a two-element vector of positive real numbers.');
end
if ~isfield(parameters,'cp_dur')
    error('CP contact time must be specified in parameters.cp_dur variable.');
elseif (~isnumeric(parameters.cp_dur))||(~isreal(parameters.cp_dur))||...
       (numel(parameters.cp_dur)~=1)||(parameters.cp_dur<=0)
    error('parameters.cp_dur must be a positive real scalar.');
end
if ~isfield(parameters,'offset')
    error('transmitter offsets must be specified in parameters.offset variable.');
elseif (~isnumeric(parameters.offset))||(~isreal(parameters.offset))||...
       (numel(parameters.offset)~=2)
    error('parameters.offset must be a two-element vector of real numbers.');
end
if ~isfield(parameters,'nblocks')
    error('the number of FSLG blocks must be specified in parameters.nblocks variable.');
elseif (~isnumeric(parameters.nblocks))||(~isreal(parameters.nblocks))||...
       (numel(parameters.nblocks)~=1)||(parameters.nblocks<1)||...
       (mod(parameters.nblocks,1)~=0)
    error('parameters.nblocks must be a positive integer.');
end
if ~isfield(parameters,'spc_dim')
    error('Fokker-Planck dimension should be specified in parameters.spc_dim variable.');
elseif (~isnumeric(parameters.spc_dim))||(~isreal(parameters.spc_dim))||...
       (numel(parameters.spc_dim)~=1)||(parameters.spc_dim<1)||...
       (mod(parameters.spc_dim,1)~=0)
    error('parameters.spc_dim should be a positive integer.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
elseif (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
       (numel(parameters.sweep)~=2)||any(parameters.sweep<=0)
    error('parameters.sweep must be a two-element vector of positive real numbers.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif (~iscell(parameters.spins))||(numel(parameters.spins)~=2)
    error('parameters.spins cell array should have exactly two elements.');
elseif any(~cellfun(@ischar,parameters.spins))
    error('elements of parameters.spins cell array must be strings.');
elseif any(~ismember(parameters.spins,spin_system.comp.isotopes))
    error('parameters.spins contains isotopes that are not present in the system.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
       (numel(parameters.npoints)~=2)||any(parameters.npoints<1)||...
       any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a two-element vector of positive integers.');
end
end

% The early bird may get the worm, but the 
% second mouse gets the cheese.
%
% Steven Wright

