% Translates Bruker acquisition parameters into a Spinach pulse
% sequence parameter structure. Syntax:
%
%             parameters=bruker2acq(acqus,point_factors)
%
% Parameters:
%
%    acqus          - Bruker acquisition parameter structure, or a cell
%                     array of structures ordered from the slowest to
%                     the directly detected dimension
%
%    point_factors  - divisors converting Bruker TD values into complex
%                     time-domain points, one per acquisition dimension
%
% Outputs:
%
%    parameters     - Spinach parameter structure containing sweep,
%                     npoints, offset, axis_units, timestep, and acq_time
%
% The input structures may use any capitalization of SW_h, SW, SFO1,
% O1, and TD field names. SW_h is preferred; SW*SFO1 is the fallback.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bruker2acq.m>

function parameters=bruker2acq(acqus,point_factors)

% Normalize the dimension list
if isstruct(acqus), acqus={acqus}; end

% Check consistency
grumble(acqus,point_factors);

% Preallocate outputs
n_dims=numel(acqus);
parameters.sweep=zeros(1,n_dims);
parameters.npoints=zeros(1,n_dims);
parameters.offset=zeros(1,n_dims);
parameters.timestep=zeros(1,n_dims);
parameters.acq_time=zeros(1,n_dims);

% Translate each acquisition dimension
for n=1:n_dims

    % Spectral width in Hz
    sweep=getpar(acqus{n},'SW_h',[]);
    if isempty(sweep)
        sweep=getpar(acqus{n},'SW',[])*getpar(acqus{n},'SFO1',[]);
    end

    % Point count and transmitter offset
    npoints=getpar(acqus{n},'TD',[])/point_factors(n);
    offset=getpar(acqus{n},'O1',[]);

    % Store Spinach acquisition settings
    parameters.sweep(n)=sweep;
    parameters.npoints(n)=npoints;
    parameters.offset(n)=offset;
    parameters.timestep(n)=1/sweep;
    parameters.acq_time(n)=(npoints-1)/sweep;

end

% Default frequency axis convention
parameters.axis_units='ppm';

end

% Case-insensitive structure field lookup
function value=getpar(source,name,default)
fields=fieldnames(source);
index=find(strcmpi(fields,name),1);
if isempty(index)
    value=default;
else
    value=source.(fields{index});
end
end

% Consistency enforcement
function grumble(acqus,point_factors)
if (~iscell(acqus))||isempty(acqus)||any(~cellfun(@isstruct,acqus))
    error('acqus must be a structure or a non-empty cell array of structures.');
end
if (~isnumeric(point_factors))||(~isreal(point_factors))||...
   (~isrow(point_factors))||(numel(point_factors)~=numel(acqus))||...
   any(~isfinite(point_factors))||any(point_factors<=0)||...
   any(mod(point_factors,1)~=0)
    error('point_factors must contain one positive integer per dimension.');
end
for n=1:numel(acqus)
    sw_h=getpar(acqus{n},'SW_h',[]);
    sw=getpar(acqus{n},'SW',[]);
    sfo1=getpar(acqus{n},'SFO1',[]);
    td=getpar(acqus{n},'TD',[]);
    o1=getpar(acqus{n},'O1',[]);
    if isempty(sw_h)&&(isempty(sw)||isempty(sfo1))
        error('each acquisition structure must contain SW_h or both SW and SFO1.');
    end
    if (~isempty(sw_h))&&((~isnumeric(sw_h))||(~isreal(sw_h))||...
       (~isscalar(sw_h))||(~isfinite(sw_h))||(sw_h<=0))
        error('SW_h values must be positive real finite scalars.');
    end
    if isempty(sw_h)&&((~isnumeric(sw))||(~isreal(sw))||(~isscalar(sw))||...
       (~isfinite(sw))||(sw<=0)||(~isnumeric(sfo1))||(~isreal(sfo1))||...
       (~isscalar(sfo1))||(~isfinite(sfo1))||(sfo1<=0))
        error('SW and SFO1 values must be positive real finite scalars.');
    end
    if isempty(td)||(~isnumeric(td))||(~isreal(td))||(~isscalar(td))||...
       (~isfinite(td))||(td<=0)||(mod(td,point_factors(n))~=0)
        error('TD values must be positive integers divisible by point_factors.');
    end
    if isempty(o1)||(~isnumeric(o1))||(~isreal(o1))||...
       (~isscalar(o1))||(~isfinite(o1))
        error('each acquisition structure must contain an O1 real finite scalar.');
    end
end
end

% The nice thing about standards is that there are so many to choose from.
%
% Andrew S. Tanenbaum


