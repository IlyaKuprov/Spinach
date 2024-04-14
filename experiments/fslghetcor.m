% Heteronuclear correlation MAS NMR experiment with frequency-switched
% Lee-Goldburg homonuclear decoupling. Further details in:
% 
%    https://doi.org/10.1103/PhysRev.140.A1261
%    https://doi.org/10.1016/0009-2614(89)87166-0
%    https://doi.org/10.1006/jmre.1996.1089
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
Hp=operator(spin_system,'L+',parameters.spins{1});
Hz=operator(spin_system,'Lz',parameters.spins{1});
Cp=operator(spin_system,'L+',parameters.spins{2});
Hp=kron(speye(parameters.spc_dim),Hp);
Hz=kron(speye(parameters.spc_dim),Hz);
Cp=kron(speye(parameters.spc_dim),Cp);
Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i; Cx=(Cp+Cp')/2;
 
% Compose the Liouvillian
L=H+1i*R+1i*K;
    
% Apply the 90-degree pulses
rho_cos=step(spin_system,L+2*pi*parameters.hi_pwr*Hx,...
             parameters.rho0,1/(4*parameters.hi_pwr));
rho_sin=step(spin_system,L+2*pi*parameters.hi_pwr*Hy,...
             parameters.rho0,1/(4*parameters.hi_pwr));

% Flip to the magic angle
rho_cos=step(spin_system,L+2*pi*parameters.hi_pwr*Hx,...
             rho_cos,acos(1/sqrt(3))/(2*pi)/parameters.hi_pwr);
rho_sin=step(spin_system,L+2*pi*parameters.hi_pwr*Hy,...
             rho_sin,acos(1/sqrt(3))/(2*pi)/parameters.hi_pwr);
 
% Get dwell times 
dwell_times=1./parameters.sweep;

% FSLG pulse evolution generators
L_FSLG_p=L-2*pi*parameters.offset(1)*Hz ...
          +2*pi*(parameters.hi_pwr/sqrt(2))*Hz;
L_FSLG_m=L-2*pi*parameters.offset(1)*Hz ...
          -2*pi*(parameters.hi_pwr/sqrt(2))*Hz;
L1_cos=L_FSLG_m+2*pi*parameters.hi_pwr*Hy;
L2_cos=L_FSLG_p-2*pi*parameters.hi_pwr*Hy;
L1_sin=L_FSLG_p+2*pi*parameters.hi_pwr*Hx;
L2_sin=L_FSLG_m-2*pi*parameters.hi_pwr*Hx;

% Backflip pulse evolution generators
L_backtoXY_cos=L-2*pi*parameters.hi_pwr*Hx;
L_backtoXY_sin=L-2*pi*parameters.hi_pwr*Hy;

% CP period evolution generator
L_c=L-2*pi*parameters.cp_pwr(1)*Hy...
     +2*pi*parameters.cp_pwr(2)*Cx;

% Preallocate fid components
fid.cos=zeros([parameters.npoints(2) ...
               parameters.npoints(1)],'like',1i);
fid.sin=zeros([parameters.npoints(2) ...
               parameters.npoints(1)],'like',1i);

% Run the bulk of the sequence
for k=1:parameters.npoints(1)

    if k==1

        % First F1 step gets rho from above
        rho_cos_ev=rho_cos; rho_sin_ev=rho_sin;

    else

        % Subsequent F1 steps keep going
        for n=1:parameters.nblocks
            rho_cos_ev=step(spin_system,L1_cos,rho_cos_ev,sqrt(2/3)/parameters.hi_pwr);
            rho_cos_ev=step(spin_system,L2_cos,rho_cos_ev,sqrt(2/3)/parameters.hi_pwr);
            rho_sin_ev=step(spin_system,L1_sin,rho_sin_ev,sqrt(2/3)/parameters.hi_pwr);
            rho_sin_ev=step(spin_system,L2_sin,rho_sin_ev,sqrt(2/3)/parameters.hi_pwr);
        end

    end
    
    % Grab the current F1 point and divert it to F2
    rho_cos_cont=rho_cos_ev; rho_sin_cont=rho_sin_ev;    
 
    % Flip the 1H back from the magic angle to the x,y-plane
    rho_cos_cont=step(spin_system,L_backtoXY_cos,rho_cos_cont,...
                      acos(1/sqrt(3))/(2*pi)/parameters.hi_pwr);
    rho_sin_cont=step(spin_system,L_backtoXY_sin,rho_sin_cont,...
                      acos(1/sqrt(3))/(2*pi)/parameters.hi_pwr);

    % Save work by stacking cos and sin parts
    rho_cont=[rho_cos_cont rho_sin_cont];

    % Run the CP contact period
    rho_cont=step(spin_system,L_c,rho_cont,parameters.cp_dur);
   
    % Decouple 1H for the acquisition period
    [L_dec,rho_cont]=decouple(spin_system,L,rho_cont,parameters.spins(1));
    
    % Run the F2 evolution
    traj=evolution(spin_system,L_dec,parameters.coil,rho_cont,...
                   dwell_times(2),parameters.npoints(2)-1,'observable');
    fid.cos(:,k)=traj(:,1); fid.sin(:,k)=traj(:,2);

end

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
if ~isfield(parameters,'hi_pwr')||(parameters.hi_pwr<=0)
    error('high RF amplitude must be specified in parameters.hi_pwr variable.');
end
if ~isfield(parameters,'cp_pwr')
    error('RF amplitude during CP must be specified in parameters.cp_pwr variable.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=2
    error('parameters.sweep array should have exactly two elements.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=2
    error('parameters.spins cell array should have exactly two elements.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('parameters.npoints array should have exactly two elements.');
end
end

% The early bird may get the worm, but the 
% second mouse gets the cheese.
%
% Steven Wright

