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
%                            within six unit roundoffs (five operations
%                            to form and recover a synchronous dwell)
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
grumble(parameters);

% Build one coherent rotor period
[prefix,~]=period_propagators(spin_system,parameters,H,R,K);
count=numel(H);

% Decompose the unitary period in an orthonormal basis
[V,T]=schur(prefix{end},'complex');
phase=angle(diag(T));
rho=V'*parameters.rho0*V;
fid=complex(zeros(1,parameters.npoints));
ticks=parameters.dwell*abs(parameters.rate)*numel(H);
nearest=round(ticks);
tolerance=(3*eps)*nearest;
if ~isfinite(ticks)||nearest<1||abs(ticks-nearest)>tolerance
    error('Acquisition dwell must be a positive integer number of rotor ticks.');
end
dwell_ticks=nearest;

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
function grumble(parameters)
validateattributes(parameters.dwell,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(parameters.npoints,{'numeric'},{'scalar','integer','positive'});
end


