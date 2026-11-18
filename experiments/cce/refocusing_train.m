% Ideal refocusing train with explicit physical preparation and detection.
% Syntax:
%
%       [answer,aux]=refocusing_train(spin_system,parameters,H,R,K)
%
% parameters.times contains elapsed seconds. refocusing is a cell row of
% ideal descriptors: physical spins, integrated Hermitian operator or
% axis/angle, and time or fraction. Alternatively train declares name (xy8
% or udd), count, axis, and spins. XY8 counts are multiples of eight; Uhrig
% fractions are sin(pi*k/(2*count+2))^2. Free-evolution preparation/detection
% conventions apply. Generalised Hilbert dynamics only. aux retains the
% compiled intervals, pulse matrices, exterior expectations, and densities.
%
% ilya.kuprov@weizmann.ac.il

function [answer,aux]=refocusing_train(spin_system,parameters,H,R,K)

% Compile experiment-owned ideal timings and propagate the retained subsystem
grumble(spin_system,parameters,H,R,K);
[answer,aux]=cce_acquire(spin_system,parameters,H,R,K);

end

% Validate the train's representation and physical control ownership
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(spin_system)||~isfield(spin_system,'formalism')||...
   ~strcmp(spin_system.formalism,'zeeman-hilb')||~isnumeric(H)||...
   ~isnumeric(R)||~isnumeric(K)||any(R(:))||any(K(:))
    error('refocusing_train:generators','Ideal trains require generalized closed Hilbert dynamics.');
end
if isfield(parameters,'sequence')
    error('refocusing_train:sequence','Use refocusing_train native descriptors, not parameters.sequence.');
end
if ~isstruct(parameters)||~isscalar(parameters)||...
   ~xor(isfield(parameters,'refocusing'),isfield(parameters,'train'))||isfield(parameters,'drive')
    error('refocusing_train:input','Supply either refocusing descriptors or one named train.');
end
end


