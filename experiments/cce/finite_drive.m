% Finite rectangular or piecewise drive in a retained cluster.
% Syntax:
%
%     [answer,aux]=finite_drive(spin_system,parameters,H,R,K)
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
% parameters.times gives increasing elapsed seconds, including during a drive.
% parameters.drive is a nonempty structure row, each element with fields:
%   spins     - unique physical observed-spin identifiers
%   axis      - 'x', 'y', 'z', or a real unit three-vector
%   amplitude - real finite amplitude row in rad/s
%   detuning  - real finite longitudinal offset in rad/s
%   duration  - matching positive duration row in seconds
% Optional time specifies the start in seconds; otherwise the start is zero.
% Alternatively fraction specifies one start as a fraction of each elapsed
% time, for a single finite pulse whose duration remains in seconds.
% Equal-length amplitude/duration rows specify consecutive waveform slices.
% A matching time row gives explicit starts, including overlapping slices.
% Optional refocusing is a cell row of ideal-pulse descriptors (spins,
% operator or axis/angle, and time or fraction). The complete drift stays
% active. Preparation, detection, output, and normalisation follow
% free_evolution. Projected branches are not supported. Pulse events
% are compiled inside this experiment.
% Mixed ideal/finite controls require generalised Hilbert dynamics.
%
% Outputs:
%
%    answer - one row per elapsed time: detected signals or observables,
%             or the column-vectorised observed density
%
%    aux    - propagation diagnostics and absolute arithmetic error bounds
%
% ilya.kuprov@weizmann.ac.il

function [answer,aux]=finite_drive(spin_system,parameters,H,R,K)

% Validate the declarative drive before compiling its finite event
grumble(spin_system,parameters,H,R,K);

% Reuse the experiment-level finite-control propagation contract
[answer,aux]=cce_acquire(spin_system,parameters,H,R,K);

end

% Validate every drive field and reject competing pulse specifications
function grumble(spin_system,parameters,H,R,K)
if ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'times')||...
   ~isnumeric(parameters.times)||~isreal(parameters.times)||~isvector(parameters.times)||...
   any(~isfinite(parameters.times))||any(parameters.times<0)||any(diff(parameters.times)<0)
    error('finite_drive:times','Supply finite nonnegative increasing elapsed times.');
end
if ~isstruct(spin_system)||~isfield(spin_system,'model')||...
   ~isstruct(parameters)||~isscalar(parameters)||~isfield(parameters,'drive')
    error('finite_drive:input','Supply a retained subsystem and drive parameters.');
end
if isfield(parameters,'drive')
    if ~isstruct(parameters.drive)||~isrow(parameters.drive)||isempty(parameters.drive)
        error('cce_pulse_events:drive','Supply a nonempty row of native drive descriptors.');
    end
    for k=1:numel(parameters.drive)
        drive=parameters.drive(k);
        if ~all(isfield(drive,{'spins','axis','amplitude','detuning','duration'}))||...
           ~isempty(setdiff(fieldnames(drive),{'spins','axis','amplitude','detuning','duration','time','fraction'}))||...
           ~isnumeric(drive.amplitude)||~isreal(drive.amplitude)||~isrow(drive.amplitude)||...
           isempty(drive.amplitude)||any(~isfinite(drive.amplitude))||...
           ~isnumeric(drive.duration)||~isreal(drive.duration)||...
           ~isequal(size(drive.duration),size(drive.amplitude))||...
           any(~isfinite(drive.duration))||any(drive.duration<=0)||...
           ~isnumeric(drive.detuning)||~isreal(drive.detuning)||...
           ~isscalar(drive.detuning)||~isfinite(drive.detuning)
            error('cce_pulse_events:drive','Supply finite row amplitudes, matching positive durations, and a real detuning.');
        end
        if ~isnumeric(drive.spins)||~isreal(drive.spins)||~isrow(drive.spins)||...
           isempty(drive.spins)||any(~isfinite(drive.spins))||any(drive.spins<1)||...
           any(drive.spins~=fix(drive.spins))||numel(unique(drive.spins))~=numel(drive.spins)
            error('cce_pulse_events:spins','Drive spins must be distinct positive physical identifiers.');
        end
        axis=drive.axis;
        if ischar(axis)&&ismember(axis,{'x','y','z'})
            axis=double('xyz'==axis);
        end
        if ~isnumeric(axis)||~isreal(axis)||~isequal(size(axis),[1 3])||...
           any(~isfinite(axis))||abs(norm(axis)-1)>3*eps
            error('cce_pulse_events:axis','A drive axis is x, y, z, or a unit three-vector.');
        end
        if isfield(drive,'time')&&(~isnumeric(drive.time)||~isreal(drive.time)||...
           ~isrow(drive.time)||(~isscalar(drive.time)&&~isequal(size(drive.time),size(drive.amplitude)))||...
           any(~isfinite(drive.time))||any(drive.time<0)||any(diff(drive.time)<0))
            error('cce_pulse_events:time','Drive starts must be finite nonnegative increasing seconds.');
        end
        if isfield(drive,'fraction')&&(isfield(drive,'time')||~isscalar(drive.amplitude)||...
           ~isnumeric(drive.fraction)||~isreal(drive.fraction)||~isscalar(drive.fraction)||...
           ~isfinite(drive.fraction)||drive.fraction<0||drive.fraction>1)
            error('cce_pulse_events:fraction','Supply one finite pulse with a start fraction in [0,1], without time.');
        end
    end
end
if isfield(parameters,'sequence')
    error('finite_drive:sequence',...
          'Use finite_drive native drive descriptors, not parameters.sequence.');
end
for k=1:numel(parameters.drive)
    ids=parameters.drive(k).spins;
    if ~isnumeric(ids)||~isreal(ids)||~isrow(ids)||isempty(ids)||...
       numel(unique(ids))~=numel(ids)||any(~ismember(ids,spin_system.model.observed_ids))
        error('finite_drive:spins',...
              'Finite drives address unique observed-spin identifiers.');
    end
end
if ~isnumeric(H)||~isnumeric(R)||~isnumeric(K)
    error('finite_drive:generators',...
          'Finite drives require conventional numeric generators.');
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


