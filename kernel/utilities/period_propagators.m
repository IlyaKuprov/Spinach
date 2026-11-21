% Builds chronological within-period prefixes of coherent MAS evolution.
% Prefix{1} is the identity; prefix{end} is the complete period.
% H, R, K, and parameters use the periodic_acquire callback convention.
% Syntax:
%
%   [prefix,dt]=period_propagators(spin_system,parameters,H,R,K)
%
% talos@spindynamics.org

function [prefix,dt]=period_propagators(spin_system,parameters,H,R,K)

% Validate coherent dynamics and a complete rotor period
grumble(spin_system,parameters,H,R,K);

% Assemble chronological ticks and rotate the combined period
count=numel(H);
order=mod(-sign(parameters.rate)*(0:count-1),count)+1;
order=circshift(order,1-parameters.start_tick);
dt=1/(abs(parameters.rate)*count);

% Compute prefixes and the one-period propagator
prefix=cell(1,count+1);
prefix{1}=eye(size(H{1}));
for n=1:count

    % Add the control belonging to the same chronological tick
    index=order(n);
    L=H{index};
    if ~isempty(parameters.tick_rf)
        L=L+parameters.tick_rf{mod(n+parameters.start_tick-2,count)+1};
    end
    prefix{n+1}=expm(-0.5i*full(L+L')*dt)*prefix{n};

end

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~strcmp(spin_system.bas.formalism,'zeeman-hilb')
    error('Periodic acquisition requires zeeman-hilb formalism.');
end
if ~isnumeric(R)||~isnumeric(K)||any(R(:)~=0)||any(K(:)~=0)
    error('Periodic acquisition requires zero relaxation and kinetics.');
end
if ~iscell(H)||~isvector(H)||isempty(H)
    error('H must contain one rotor period of Hamiltonians.');
end
shape=size(H{1});
if numel(shape)~=2||shape(1)~=shape(2)||...
   any(~cellfun(@(x) isnumeric(x)&&isequal(size(x),shape)&&...
                    norm(x-x','fro')<=spin_system.tols.liouv_zero*numel(x),H))
    error('The rotor Hamiltonians must be equally sized Hermitian matrices.');
end
if ~iscell(parameters.tick_rf)||(~isempty(parameters.tick_rf)&&...
   (~isvector(parameters.tick_rf)||numel(parameters.tick_rf)~=numel(H)||...
    any(~cellfun(@(x) isnumeric(x)&&isequal(size(x),shape)&&ishermitian(x),parameters.tick_rf))))
    error('Supply either no RF or exactly one period of Hermitian RF matrices.');
end
if ~isequal(size(parameters.rho0),size(parameters.coil),shape)
    error('Initial state and receiver must match the Hamiltonian dimensions.');
end
validateattributes(parameters.rate,{'numeric'},{'scalar','real','finite','nonzero'});
validateattributes(parameters.start_tick,{'numeric'},{'scalar','integer','positive','<=',numel(H)});
end


