% Phase-sensitive heteronuclear NOESY pulse sequence. Syntax:
%
%            fid=hoesy(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.sweep     two sweep widths, Hz
%
%     parameters.npoints   number of FID points for both
%                          dimensions
%
%     parameters.spins     nuclei on which the sequence runs,
%                          e.g. {'15N','13C'}
%
%     parameters.tmix      mixing time, seconds
%
%     parameters.needs     should be set to {'rho_eq'}, this
%                          sequence needs the thermal equili-
%                          brium state
%
%     H - Hamiltonian matrix, received from context function
%
%     R - relaxation superoperator, received from context function
%
%     K - kinetics superoperator, received from context function
%
% Outputs:
%
%     fid.cos,fid.sin - two components of the FID for F1 hyper-
%                       complex processing
%
% Zak El-Machachi, Ilya Kuprov
%
% <https://spindynamics.org/wiki/index.php?title=hoesy.m>

function fid=hoesy(spin_system,parameters,H,R,K)

% Consistency check
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Coherent evolution timestep
timestep=1./parameters.sweep;

% Detection state
coil=state(spin_system,'L+',parameters.spins{2},'cheap');

% Pulse operators
Hp=operator(spin_system,'L+',parameters.spins{1});
Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i;
Cp=operator(spin_system,'L+',parameters.spins{2});
Cy=(Cp-Cp')/2i;

% First pulse on F1
rho=step(spin_system,Hx,parameters.rho0,pi/2);

% First half of F1 evolution
rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2,...
                    parameters.npoints(1)-1,'trajectory');

% F1 dimension decoupling
for n=1:numel(parameters.decouple_f1)
    Lp=operator(spin_system,'L+',parameters.decouple_f1{n});
    rho_stack=step(spin_system,(Lp+Lp')/2,rho_stack,pi);
end

% Second half of F1 evolution
rho_stack=evolution(spin_system,L,[],rho_stack,timestep(1)/2,...
                    parameters.npoints(1)-1,'refocus');

% Second pulse on F1
rho_stack_cos_p=step(spin_system,Hx,rho_stack,+pi/2);
rho_stack_sin_p=step(spin_system,Hy,rho_stack,+pi/2);
rho_stack_cos_m=step(spin_system,Hx,rho_stack,-pi/2);
rho_stack_sin_m=step(spin_system,Hy,rho_stack,-pi/2);

% Homospoil
rho_stack_cos_p=homospoil(spin_system,rho_stack_cos_p,'destroy');
rho_stack_sin_p=homospoil(spin_system,rho_stack_sin_p,'destroy');
rho_stack_cos_m=homospoil(spin_system,rho_stack_cos_m,'destroy');
rho_stack_sin_m=homospoil(spin_system,rho_stack_sin_m,'destroy');

% Mixing time
rho_stack_cos_p=evolution(spin_system,1i*R+1i*K,[],rho_stack_cos_p,parameters.tmix,1,'final');
rho_stack_sin_p=evolution(spin_system,1i*R+1i*K,[],rho_stack_sin_p,parameters.tmix,1,'final');
rho_stack_cos_m=evolution(spin_system,1i*R+1i*K,[],rho_stack_cos_m,parameters.tmix,1,'final');
rho_stack_sin_m=evolution(spin_system,1i*R+1i*K,[],rho_stack_sin_m,parameters.tmix,1,'final');

% Homospoil
rho_stack_cos_p=homospoil(spin_system,rho_stack_cos_p,'destroy');
rho_stack_sin_p=homospoil(spin_system,rho_stack_sin_p,'destroy');
rho_stack_cos_m=homospoil(spin_system,rho_stack_cos_m,'destroy');
rho_stack_sin_m=homospoil(spin_system,rho_stack_sin_m,'destroy');

% Pulse on F2
rho_stack_cos_p=step(spin_system,Cy,rho_stack_cos_p,pi/2);
rho_stack_sin_p=step(spin_system,Cy,rho_stack_sin_p,pi/2);
rho_stack_cos_m=step(spin_system,Cy,rho_stack_cos_m,pi/2);
rho_stack_sin_m=step(spin_system,Cy,rho_stack_sin_m,pi/2);

% Axial peak elimination in F2
rho_stack_cos=rho_stack_cos_p-rho_stack_cos_m;
rho_stack_sin=rho_stack_sin_p-rho_stack_sin_m;

% Decouple protons
[L,rho_stack_cos]=decouple(spin_system,L,rho_stack_cos,parameters.spins(1));
[L,rho_stack_sin]=decouple(spin_system,L,rho_stack_sin,parameters.spins(1));

% F2 evolution
fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep(2),parameters.npoints(2)-1,'observable');
fid.sin=evolution(spin_system,L,coil,rho_stack_sin,timestep(2),parameters.npoints(2)-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only available for sphten-liouv and zeeman-liouv formalisms.');
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
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
elseif numel(parameters.spins)~=2
    error('parameters.spins cell array should have exactly one element.');
end
if ~isfield(parameters,'npoints')
    error('number of points should be specified in parameters.npoints variable.');
elseif numel(parameters.npoints)~=2
    error('parameters.npoints array should have exactly two elements.');
end
if ~isfield(parameters,'tmix')
    error('mixing time should be specified in parameters.tmix variable.');
elseif numel(parameters.tmix)~=1
    error('parameters.tmix array should have exactly one element.');
end
end

% Divide by cucumber error. Please 
% reinstall Universe and reboot.
%
% Terry Pratchett

