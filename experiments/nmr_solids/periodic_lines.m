% Resolves the per-crystal line spectrum of a coherent MAS period.
% Set sum_up=false in singlerot: line lists do not share a frequency
% axis. Inputs are as in periodic_acquire, without dwell or npoints.
% Each row contains frequency (Hz), complex amplitude, and sideband.
% Syntax: lines=periodic_lines(spin_system,parameters,H,R,K)
%
% talos@spindynamics.org

function lines=periodic_lines(spin_system,parameters,H,R,K)

% Check that orientations will remain separate
grumble(parameters);

% Resolve the period spectrum without constructing a finite FID
[prefix,~]=period_propagators(spin_system,parameters,H,R,K);
[freq,amp,band]=period_lines(prefix,abs(parameters.rate),parameters.rho0,parameters.coil);
lines=[freq amp band];

end

% Consistency enforcement
function grumble(parameters)
if ~isfield(parameters,'sum_up')||~isequal(parameters.sum_up,false)
    error('Periodic line acquisition requires parameters.sum_up=false.');
end
end


