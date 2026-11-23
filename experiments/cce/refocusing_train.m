% Ideal refocusing train with explicit physical preparation and detection.
% Syntax:
%
%       [answer,aux]=refocusing_train(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    spin_system - native retained subsystem from cce_subsystem
%
%    parameters  - experiment structure with the fields described below
%
%    H,R,K       - angular-frequency generators in the retained basis;
%                  L=H+1i*R+1i*K for conventional matrix generators
%
% parameters.times contains elapsed seconds. refocusing is a cell row of
% ideal descriptors: physical spins, integrated Hermitian operator or
% axis/angle, and time or fraction. Alternatively train declares name (xy8
% or udd), count, axis, and spins. XY8 counts are multiples of eight; Uhrig
% fractions are sin(pi*k/(2*count+2))^2. Free-evolution preparation/detection
% conventions apply. Generalised Hilbert dynamics only. aux retains the
% compiled intervals, pulse matrices, exterior expectations, and densities.
%
% Outputs:
%
%    answer - one row per elapsed time: detected signals or observables,
%             or the column-vectorised observed density
%
%    aux    - arithmetic error bounds, compiled segments, and densities
%
% ilya.kuprov@weizmann.ac.il

function [answer,aux]=refocusing_train(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compile ideal pulse timings and propagate the retained subsystem
[answer,aux]=cce_acquire(spin_system,parameters,H,R,K);

end

% Validate the train's representation and physical control ownership
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'times')||...
   ~isnumeric(parameters.times)||~isreal(parameters.times)||~isvector(parameters.times)||...
   any(~isfinite(parameters.times))||any(parameters.times<0)||any(diff(parameters.times)<0)
    error('refocusing_train:times','Supply finite nonnegative increasing elapsed times.');
end
if ~isstruct(spin_system)||~isfield(spin_system,'bas')||...
   ~strcmp(spin_system.bas.formalism,'zeeman-hilb')||~isnumeric(H)||...
   ~isnumeric(R)||~isnumeric(K)||any(R(:))||any(K(:))
    error('refocusing_train:generators',...
          'Ideal trains require generalized closed Hilbert dynamics.');
end
if isfield(parameters,'sequence')
    error('refocusing_train:sequence',...
          'Use refocusing_train native descriptors, not parameters.sequence.');
end
if ~isstruct(parameters)||~isscalar(parameters)||...
   ~xor(isfield(parameters,'refocusing'),isfield(parameters,'train'))||isfield(parameters,'drive')
    error('refocusing_train:input',...
          'Supply either refocusing descriptors or one named train.');
end
if isfield(parameters,'refocusing')
    if ~iscell(parameters.refocusing)||~isrow(parameters.refocusing)
        error('cce_pulse_events:refocusing','Supply a cell row of ideal-pulse descriptors.');
    end
    for k=1:numel(parameters.refocusing)
        pulse=parameters.refocusing{k};
        if ~isstruct(pulse)||~isscalar(pulse)||~isfield(pulse,'spins')||...
           ~xor(isfield(pulse,'operator'),all(isfield(pulse,{'axis','angle'})))||...
           ~xor(isfield(pulse,'fraction'),isfield(pulse,'time'))||...
           ~isempty(setdiff(fieldnames(pulse),{'spins','operator','axis','angle','fraction','time'}))
            error('cce_pulse_events:refocusing','Supply spins, integrated operator or axis/angle, and time or fraction.');
        end
    end
end
if isfield(parameters,'train')
    train=parameters.train;
    if ~isstruct(train)||~isscalar(train)||...
       ~all(isfield(train,{'name','count','axis','spins'}))||...
       ~isempty(setdiff(fieldnames(train),{'name','count','axis','spins'}))||...
       ~ismember(train.name,{'xy8','udd'})||~ischar(train.axis)||...
       ~ismember(train.axis,{'x','y','z'})||~isnumeric(train.count)||...
       ~isreal(train.count)||~isscalar(train.count)||~isfinite(train.count)||...
       train.count<1||train.count~=fix(train.count)||...
       (strcmp(train.name,'xy8')&&mod(train.count,8)~=0)
        error('cce_pulse_events:train','Supply an XY8 multiple of eight or a positive Uhrig pulse count.');
    end
end
end


