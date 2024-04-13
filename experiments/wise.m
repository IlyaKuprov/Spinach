% WISE (WIdeline SEparation) is a powder MAS heteronuclear correlation 
% experiment. In the common 1H-13C implementation, molecular dynamics
% information is contained in 1H line shapes that are separated in the
% second dimension by 13C chemical shifts. Further information in:
% 
%                 https://doi.org/10.1021/ma00038a037
%
% Syntax:
%
%                fid=wise(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.spins        - working spins, e.g. {'1H',13C'}
%
%    parameters.hi_pwr       - amplitude of high power pulses 
%                              on the high-gamma channel, Hz
%
%    parameters.cp_pwr       - amplitude of CP pulse on each 
%                              channel during the CP contact
%                              time, Hz
%
%    parameters.cp_dur       - CP contact time duration, s
%
%    parameters.rho0         - initial state
%
%    parameters.coil         - detection state
%
%    parameters.sweep        - sweep width, Hz for F1, F2
%
%    parameters.npoints      - number of points in F1, F2
% 
%    H  - Hamiltonian matrix, received from context function
%
%    R  - relaxation superoperator, received from context function
%
%    K  - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid.sin, fid.cos        - sine and cosine components
%                              of the States quadrature
%
% guinevere.mathies@uni-konstanz.de
%
% <https://spindynamics.org/wiki/index.php?title=wise.m>

function fid=wise(spin_system,parameters,H,R,K)
 
% Consistency enforcement
grumble(spin_system,parameters,H,R,K);

% Wipe the state of 13C (pre-saturation)
[~,parameters.rho0]=decouple(spin_system,[],parameters.rho0,{'13C'});
 
% Build 1H and 13C control operators
Hp=operator(spin_system,'L+',parameters.spins{1});
Cp=operator(spin_system,'L+',parameters.spins{2});
Hp=kron(speye(parameters.spc_dim),Hp);
Cp=kron(speye(parameters.spc_dim),Cp);
Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i; Cx=(Cp+Cp')/2;
 
% Compose the Liouvillian
L=H+1i*R+1i*K;
    
% High-power 90-degree pulses on 1H along X (cos) and Y (sin)
L_hp_cos=L+2*pi*parameters.hi_pwr*Hx; L_hp_sin=L+2*pi*parameters.hi_pwr*Hy;
rho_cos=step(spin_system,L_hp_cos,parameters.rho0,1/(4*parameters.hi_pwr));
rho_sin=step(spin_system,L_hp_sin,parameters.rho0,1/(4*parameters.hi_pwr));
 
% Get dwell times
dw=1./parameters.sweep;
 
% Run the F1 evolution
rho_stack_cos=evolution(spin_system,L,[],rho_cos,dw(1),...
                        parameters.npoints(1)-1,'trajectory');
rho_stack_sin=evolution(spin_system,L,[],rho_sin,dw(1),...
                        parameters.npoints(1)-1,'trajectory');

% CP contact time evolution generator (-Y on 1H, +X on 13C)
L_cp=L-2*pi*parameters.cp_pwr(1)*Hy+2*pi*parameters.cp_pwr(2)*Cx;
   
% Run CP contact time evolution
rho_stack_cos=evolution(spin_system,L_cp,[],rho_stack_cos,...
                        parameters.cp_dur,1,'final');
rho_stack_sin=evolution(spin_system,L_cp,[],rho_stack_sin,...
                        parameters.cp_dur,1,'final');
 
% Wipe and decouple protons for acquisition
[L_dec,rho_stack_cos]=decouple(spin_system,L,rho_stack_cos,parameters.spins(1));
[~,rho_stack_sin]=decouple(spin_system,[],rho_stack_sin,parameters.spins(1));

% Run the F2 evolution
fid.cos=evolution(spin_system,L_dec,parameters.coil,rho_stack_cos,...
                  dw(2),parameters.npoints(2)-1,'observable');
fid.sin=evolution(spin_system,L_dec,parameters.coil,rho_stack_sin,...
                  dw(2),parameters.npoints(2)-1,'observable');

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

% [...] but all that pales next to the fundamental source of 
% opposition: our resentment of being pushed around by people
% who tell us it's for our own good.
%
% Peter Wood

