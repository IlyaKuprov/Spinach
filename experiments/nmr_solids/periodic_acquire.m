% Samples coherent periodic MAS evolution at integer rotor ticks.
% Use as a zeeman-hilb singlerot callback. Hamiltonians and optional
% RF contributions are piecewise constant, in chronological order
% starting at the supplied rotor phase. A start tick rotates both
% lists together; it does not rotate the state or the receiver.
%
% Syntax:
%
%           fid=periodic_acquire(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.rate       - non-zero signed rotor rate in Hz
%    parameters.tick_rf    - empty cell or one square RF matrix per
%                            chronological tick, in rad/s
%    parameters.start_tick - one-based index of the first tick
%    parameters.dwell      - acquisition dwell in seconds, an
%                            integer multiple of the rotor tick
%    parameters.npoints    - positive number of acquired samples
%    parameters.rho0       - initial Hilbert-space density matrix
%    parameters.coil       - Hilbert-space detection matrix
%    H                     - cell vector of Hamiltonians in rad/s
%                            in Fourier grid order from singlerot
%    R,K                   - zero relaxation and kinetics matrices
%
% Outputs:
%
%    fid                   - complex row, starting at time zero
%
% The generator at each Fourier node is held for one tick, starting
% at that node, with decreasing phase for positive rates. This is
% a left-endpoint rule; refine the rotor grid for continuous MAS.
%
% talos@spindynamics.org

function fid=periodic_acquire(spin_system,parameters,H,R,K)

% Check consistency
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

% Decompose the unitary period in an orthonormal basis
[V,T]=schur(prefix{end},'complex');
phase=angle(diag(T));
rho=V'*parameters.rho0*V;
fid=complex(zeros(1,parameters.npoints));
dwell_ticks=round(parameters.dwell/dt);

% Combine the integer period power with the residual prefix
for n=1:numel(fid)
    ticks=(n-1)*dwell_ticks;
    periods=floor(ticks/count);
    U=prefix{mod(ticks,count)+1}*V;
    phase_factor=exp(1i*periods*phase);
    evolved=U*(rho.*(phase_factor*phase_factor'))*U';
    fid(n)=trace(parameters.coil'*evolved);
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
validateattributes(parameters.dwell,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(parameters.npoints,{'numeric'},{'scalar','integer','positive'});
ticks=parameters.dwell*abs(parameters.rate)*numel(H);
if ticks<1||abs(ticks-round(ticks))>eps(ticks)
    error('Acquisition dwell must be a positive integer number of rotor ticks.');
end
end


