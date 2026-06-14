% Sensitivity-improved echo/antiecho gradient-selected HSQC pulse
% sequence, based on the Bruker hsqcetgpsi pulse program and the
% standard HSQC sequence from:
%
%           https://doi.org/10.1016/0009-2614(80)80041-8
%           https://doi.org/10.1002/cmr.a.10095
%           https://doi.org/10.1016/0022-2364(91)90036-S
%           https://doi.org/10.1021/ja00052a088
%           https://doi.org/10.1007/BF00175254
%
% The gradient selection is represented analytically by coherence order
% selection statements
%
% Syntax:
%
%              fid=hsqcetgpsi(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.sweep              [F1 F2] sweep widths, Hz
%
%     parameters.npoints            [F1 F2] numbers of points
%
%     parameters.spins              {F1 F2} nuclei (e.g. '13C','1H')
%
%     parameters.decouple_f2        nuclei to decouple in F2, e.g.
%                                   {'15N','13C'}
%
%     parameters.decouple_f1        nuclei that receive midpoint
%                                   180-degree refocusing pulses in
%                                   F1, e.g. {'1H','13C'}
%
%     parameters.J                  working scalar coupling, Hz
%
%     parameters.trim_angle         proton trim pulse angle, rad
%
%     parameters.si_time            sensitivity improvement delay, s
%
%     H  - Hamiltonian matrix, received from context function
%
%     R  - relaxation superoperator, received from context function
%
%     K  - kinetics superoperator, received from context function
%
% Outputs:
%
%     fid.pos,fid.neg -  echo and antiecho components of the
%                        signal.
%
% Note: natural abundance simulations should make use of the isotope
%       dilution functionality. See dilute.m function.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=hsqcetgpsi.m>

function fid=hsqcetgpsi(spin_system,parameters,H,R,K)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Coherent evolution timesteps
timestep=1./parameters.sweep;

% J-coupling evolution time
delta=abs(1/(2*parameters.J));

% Initial condition
rho=state(spin_system,'Lz',parameters.spins{2},'cheap');

% Detection state
coil=state(spin_system,'L+',parameters.spins{2},'cheap');

% Pulse operators
Cx=operator(spin_system,'Lx',parameters.spins{1});
Cy=operator(spin_system,'Ly',parameters.spins{1});
Hx=operator(spin_system,'Lx',parameters.spins{2});
Hy=operator(spin_system,'Ly',parameters.spins{2});

% Initial proton pulse
rho=step(spin_system,Hx,rho,pi/2);

% First INEPT evolution
rho=evolution(spin_system,L,[],rho,delta/2,1,'final');

% Refocusing pulses
rho=step(spin_system,Hx+Cx,rho,pi);

% Second INEPT evolution
rho=evolution(spin_system,L,[],rho,delta/2,1,'final');

% Proton trim pulse
rho=step(spin_system,Hx,rho,parameters.trim_angle);

% Transfer pulses
rho=step(spin_system,Hy,rho,pi/2);
rho=step(spin_system,Cx,rho,pi/2)-step(spin_system,Cx,rho,-pi/2);

% First half of F1 evolution
rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2,...
                    parameters.npoints(1)-1,'trajectory');

% F1 dimension refocusing pulses
for n=1:numel(parameters.decouple_f1)
    Lp=operator(spin_system,'L+',parameters.decouple_f1{n});
    rho_stack=step(spin_system,(Lp+Lp')/2,rho_stack,pi);
end

% Second half of F1 evolution
rho_stack=evolution(spin_system,L,[],rho_stack,timestep(1)/2,...
                    parameters.npoints(1)-1,'refocus');

% First analytical gradient selection
rho_stack_pos=coherence(spin_system,rho_stack,{{parameters.spins{2},0},...
                                               {parameters.spins{1},+1}});
rho_stack_neg=coherence(spin_system,rho_stack,{{parameters.spins{2},0},...
                                               {parameters.spins{1},-1}});

% Carbon inversion pulse after the first gradient
rho_stack_pos=step(spin_system,Cx,rho_stack_pos,pi);
rho_stack_neg=step(spin_system,Cx,rho_stack_neg,pi);

% First sensitivity improvement pulse pair
rho_stack_pos=step(spin_system,Hx,rho_stack_pos,pi/2);
rho_stack_pos=step(spin_system,Cx,rho_stack_pos,pi/2);
rho_stack_neg=step(spin_system,Hx,rho_stack_neg,pi/2);
rho_stack_neg=step(spin_system,Cx,rho_stack_neg,pi/2);

% First sensitivity improvement evolution period
rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,...
                        parameters.si_time,1,'final');
rho_stack_neg=evolution(spin_system,L,[],rho_stack_neg,...
                        parameters.si_time,1,'final');

% Sensitivity improvement refocusing pulses
rho_stack_pos=step(spin_system,Hx+Cx,rho_stack_pos,pi);
rho_stack_neg=step(spin_system,Hx+Cx,rho_stack_neg,pi);

% Second sensitivity improvement evolution period
rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,...
                        parameters.si_time,1,'final');
rho_stack_neg=evolution(spin_system,L,[],rho_stack_neg,...
                        parameters.si_time,1,'final');

% Second sensitivity improvement pulse pair
rho_stack_pos=step(spin_system,Hy,rho_stack_pos,pi/2);
rho_stack_pos=step(spin_system,Cy,rho_stack_pos,pi/2);
rho_stack_neg=step(spin_system,Hy,rho_stack_neg,pi/2);
rho_stack_neg=step(spin_system,Cy,rho_stack_neg,pi/2);

% First back-transfer evolution period
rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,delta/2,1,'final');
rho_stack_neg=evolution(spin_system,L,[],rho_stack_neg,delta/2,1,'final');

% Back-transfer refocusing pulses
rho_stack_pos=step(spin_system,Hx+Cx,rho_stack_pos,pi);
rho_stack_neg=step(spin_system,Hx+Cx,rho_stack_neg,pi);

% Second back-transfer evolution period
rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,delta/2,1,'final');
rho_stack_neg=evolution(spin_system,L,[],rho_stack_neg,delta/2,1,'final');

% Final proton pulse
rho_stack_pos=step(spin_system,Hx,rho_stack_pos,pi/2);
rho_stack_neg=step(spin_system,Hx,rho_stack_neg,pi/2);

% Proton echo pulse before acquisition
rho_stack_pos=step(spin_system,Hx,rho_stack_pos,pi);
rho_stack_neg=step(spin_system,Hx,rho_stack_neg,pi);

% Second analytical gradient selection
rho_stack_pos=coherence(spin_system,rho_stack_pos,{{parameters.spins{1},0},...
                                                   {parameters.spins{2},+1}});
rho_stack_neg=coherence(spin_system,rho_stack_neg,{{parameters.spins{1},0},...
                                                   {parameters.spins{2},+1}});

% Decoupling in F2
[L,rho_stack_pos]=decouple(spin_system,L,rho_stack_pos,parameters.decouple_f2);
[L,rho_stack_neg]=decouple(spin_system,L,rho_stack_neg,parameters.decouple_f2);

% Detection
fid.pos=evolution(spin_system,L,coil,rho_stack_pos,...
                  timestep(2),parameters.npoints(2)-1,'observable');
fid.neg=evolution(spin_system,L,coil,rho_stack_neg,...
                  timestep(2),parameters.npoints(2)-1,'observable');

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
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
elseif numel(parameters.sweep)~=2
    error('parameters.sweep array should have exactly two elements.');
elseif (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
       any(~isfinite(parameters.sweep))||any(parameters.sweep<=0)
    error('parameters.sweep must contain two positive real numbers.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=2
    error('parameters.spins cell array should have exactly two elements.');
elseif (~iscell(parameters.spins))||(~ischar(parameters.spins{1}))||...
       (~ischar(parameters.spins{2}))
    error('parameters.spins must be a two-element cell array of character strings.');
elseif strcmp(parameters.spins{1},parameters.spins{2})
    error('parameters.spins must specify two different isotopes.');
elseif any(~ismember(parameters.spins,spin_system.comp.isotopes))
    error('parameters.spins contains isotopes that are not present in the system.');
end
if ~isfield(parameters,'decouple_f2')
    error('decoupling channel list should be specified in parameters.decouple_f2 variable.');
elseif ~iscell(parameters.decouple_f2)
    error('parameters.decouple_f2 must be a cell array.');
elseif any(~ismember(parameters.decouple_f2,spin_system.comp.isotopes))
    error('parameters.decouple_f2 contains isotopes that are not present in the system.');
end
if ~isfield(parameters,'decouple_f1')
    error('decoupling channel list should be specified in parameters.decouple_f1 variable.');
elseif ~iscell(parameters.decouple_f1)
    error('parameters.decouple_f1 must be a cell array.');
elseif any(~ismember(parameters.decouple_f1,spin_system.comp.isotopes))
    error('parameters.decouple_f1 contains isotopes that are not present in the system.');
elseif ismember(parameters.spins{1},parameters.decouple_f1)
    error('parameters.decouple_f1 must not contain the active indirect isotope.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('parameters.npoints array should have exactly two elements.');
elseif (~isnumeric(parameters.npoints))||(~isreal(parameters.npoints))||...
       any(parameters.npoints<1)||any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must contain two positive integers.');
end
if ~isfield(parameters,'J')
    error('scalar coupling should be specified in parameters.J variable.');
elseif numel(parameters.J)~=1
    error('parameters.J array should have exactly one element.');
elseif (~isnumeric(parameters.J))||(~isreal(parameters.J))||...
       (~isfinite(parameters.J))||(parameters.J==0)
    error('parameters.J must be a non-zero real scalar.');
end
if ~isfield(parameters,'trim_angle')
    error('trim pulse angle should be specified in parameters.trim_angle variable.');
elseif numel(parameters.trim_angle)~=1
    error('parameters.trim_angle array should have exactly one element.');
elseif (~isnumeric(parameters.trim_angle))||(~isreal(parameters.trim_angle))||...
       (~isfinite(parameters.trim_angle))
    error('parameters.trim_angle must be a finite real scalar.');
end
if ~isfield(parameters,'si_time')
    error('sensitivity improvement delay should be specified in parameters.si_time variable.');
elseif numel(parameters.si_time)~=1
    error('parameters.si_time array should have exactly one element.');
elseif (~isnumeric(parameters.si_time))||(~isreal(parameters.si_time))||...
       (~isfinite(parameters.si_time))||(parameters.si_time<=0)
    error('parameters.si_time must be a positive real scalar.');
end
end


