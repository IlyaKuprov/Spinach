% Compile native experiment declarations into internal pulse events.
% Syntax:
%
%                    events=cce_pulse_events(parameters)
%
% Internal service for native acquisition. Each drive structure contains
% spins, axis, amplitude (rad/s), detuning (rad/s), and duration (seconds).
% Equal-length row amplitudes/durations describe consecutive rectangular
% waveform slices. Optional time is the pulse start in seconds, otherwise
% zero; a matching time row specifies explicit slice starts. A single pulse
% may instead declare fraction in [0,1]. A drive structure row adds independently
% phased, timed, and detuned pulses. refocusing is a cell row of physical ideal-pulse descriptors: spins,
% either integrated Hermitian operator or axis/angle, and time or fraction.
% train contains name (xy8 or udd), count, axis, and spins. No raw event
% array, omitted operator angle, or legacy parameter spelling is accepted.
%
% Parameters:
%
%   parameters - scalar native experiment structure containing drive,
%                refocusing, or train declarations in the units above
%
% Outputs:
%
%   events - structure row of finite and ideal pulse events with physical
%            targets, absolute/fractional placements, and pulse parameters
%
% ilya.kuprov@weizmann.ac.il

function events=cce_pulse_events(parameters)

% Validate declarations before generating an internal event array
grumble(parameters);
count=0;
finite_count=0;
if isfield(parameters,'drive')
    finite_count=sum(arrayfun(@(d)numel(d.amplitude),parameters.drive));
end
if isfield(parameters,'refocusing')
    count=numel(parameters.refocusing);
end
if isfield(parameters,'train')
    count=count+parameters.train.count;
end
events=repmat(struct('type','','target',[],'placement','','at',[],'duration',[],...
                     'axis',[],'amplitude',[],'detuning',[],'operator',[],'angle',[]),...
              1,finite_count+count);

% Expand the declared finite waveform without altering its sample values
if isfield(parameters,'drive')
    index=0;
    for n=1:numel(parameters.drive)
        drive=parameters.drive(n);
        start=0;
        placement='absolute';
        if isfield(drive,'time')
            start=drive.time(1);
        end
        if isfield(drive,'fraction')
            start=drive.fraction;
            placement='fraction';
        end
        for k=1:numel(drive.amplitude)
            if isfield(drive,'time')&&numel(drive.time)>1
                start=drive.time(k);
            end
            index=index+1;
            events(index)=struct('type','finite','target',drive.spins,...
                                 'placement',placement,'at',start,...
                                 'duration',drive.duration(k),'axis',drive.axis,...
                                 'amplitude',drive.amplitude(k),'detuning',drive.detuning,...
                                 'operator',[],'angle',[]);
            start=start+drive.duration(k);
        end
    end
end

% Resolve literal XY8 phases and Uhrig fractions inside the experiment
pulses={};
if isfield(parameters,'refocusing')
    pulses=parameters.refocusing;
end
if isfield(parameters,'train')
    train=parameters.train;
    pulses=cell(1,train.count);
    if strcmp(train.name,'xy8')
        phases=repmat('xyxyyxyx',1,train.count/8);
        fractions=((1:train.count)-0.5)/train.count;
    else
        phases=repmat(train.axis,1,train.count);
        fractions=sin(pi*(1:train.count)/(2*train.count+2)).^2;
    end
    for k=1:train.count
        pulses{k}=struct('spins',train.spins,'axis',phases(k),...
                         'angle',pi,'fraction',fractions(k));
    end
end

% Retain input order for coincident ideal pulses
for k=1:numel(pulses)
    pulse=pulses{k};
    placement='absolute';
    if isfield(pulse,'fraction')
        at=pulse.fraction;
        placement='fraction';
    else
        at=pulse.time;
    end
    operator=[];
    axis=[];
    angle=1;
    if isfield(pulse,'operator'), operator=pulse.operator;
    else
        axis=pulse.axis;
        angle=pulse.angle;
    end
    events(finite_count+k)=struct('type','ideal','target',pulse.spins,...
                                 'placement',placement,'at',at,'duration',[],'axis',axis,...
                                 'amplitude',[],'detuning',[],'operator',operator,'angle',angle);
end

end

% Require explicit physical declarations and reject competing old inputs
function grumble(parameters)
if ~isstruct(parameters)||~isscalar(parameters)||...
   isfield(parameters,'sequence')
    error('cce_pulse_events:input','Use finite_drive, refocusing_train, or hahn_echo native descriptors, not parameters.sequence.');
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


