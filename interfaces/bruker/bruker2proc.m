% Translates Bruker processing parameters into fields used by Spinach
% examples for apodisation and Fourier transformation. Syntax:
%
%            parameters=bruker2proc(procs,parameters)
%
% Parameters:
%
%    procs       - Bruker processing parameter structure, or a cell
%                  array of structures ordered from the slowest to
%                  the directly detected dimension
%
%    parameters  - Spinach pulse sequence parameter structure produced
%                  by bruker2acq()
%
% Outputs:
%
%    parameters  - input structure augmented with zerofill and
%                  apodisation fields; offset is updated from OFFSET,
%                  SF, and SW_p when those fields are present
%
% Bruker WDW values 0 (no window) and 1 (exponential) are supported.
% Exponential LB values in Hz are converted into the dimensionless
% coefficient expected by Spinach apodisation().
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bruker2proc.m>

function parameters=bruker2proc(procs,parameters)

% Normalize the dimension list
if isstruct(procs), procs={procs}; end

% Check consistency
grumble(procs,parameters);

% Preallocate outputs
n_dims=numel(procs);
parameters.zerofill=zeros(1,n_dims);
parameters.apodisation=cell(1,n_dims);

% Translate each processing dimension
for n=1:n_dims

    % Fourier transform size
    parameters.zerofill(n)=getpar(procs{n},'SI',[]);

    % Processed spectral axis centre
    axis_edge=getpar(procs{n},'OFFSET',[]);
    obs_freq=getpar(procs{n},'SF',[]);
    proc_sweep=getpar(procs{n},'SW_p',[]);
    if (~isempty(axis_edge))&&(~isempty(obs_freq))&&(~isempty(proc_sweep))
        parameters.offset(n)=axis_edge*obs_freq-proc_sweep/2;
    end

    % Apodisation specification
    window=getpar(procs{n},'WDW',0);
    switch window
        case 0
            parameters.apodisation{n}={'none'};
        case 1
            line_broadening=getpar(procs{n},'LB',0);
            parameters.apodisation{n}={'exp',pi*line_broadening*...
                                            parameters.acq_time(n)};
        otherwise
            error('only Bruker WDW values 0 and 1 are currently supported.');
    end

end

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
function grumble(procs,parameters)
if (~iscell(procs))||isempty(procs)||any(~cellfun(@isstruct,procs))
    error('procs must be a structure or a non-empty cell array of structures.');
end
if (~isstruct(parameters))||(~isfield(parameters,'sweep'))||...
   (~isfield(parameters,'npoints'))||(~isfield(parameters,'offset'))||...
   (~isfield(parameters,'acq_time'))
    error('parameters must be a structure produced by bruker2acq().');
end
if (numel(parameters.sweep)~=numel(procs))||...
   (numel(parameters.npoints)~=numel(procs))||...
   (numel(parameters.offset)~=numel(procs))||...
   (numel(parameters.acq_time)~=numel(procs))
    error('the number of processing dimensions must match the acquisition dimensions.');
end
for n=1:numel(procs)
    si=getpar(procs{n},'SI',[]);
    wdw=getpar(procs{n},'WDW',0);
    lb=getpar(procs{n},'LB',0);
    axis_edge=getpar(procs{n},'OFFSET',[]);
    obs_freq=getpar(procs{n},'SF',[]);
    proc_sweep=getpar(procs{n},'SW_p',[]);
    if isempty(si)||(~isnumeric(si))||(~isreal(si))||(~isscalar(si))||...
       (~isfinite(si))||(si<parameters.npoints(n))||(mod(si,1)~=0)
        error('SI values must be positive integers not smaller than npoints.');
    end
    if (~isnumeric(wdw))||(~isreal(wdw))||(~isscalar(wdw))||...
       (~isfinite(wdw))||(mod(wdw,1)~=0)
        error('WDW values must be real finite integer scalars.');
    end
    if (~isnumeric(lb))||(~isreal(lb))||(~isscalar(lb))||...
       (~isfinite(lb))||(lb<0)
        error('LB values must be non-negative real finite scalars.');
    end
    axis_fields={axis_edge,obs_freq,proc_sweep};
    if any(cellfun(@isempty,axis_fields))&&(~all(cellfun(@isempty,axis_fields)))
        error('OFFSET, SF, and SW_p must either all be present or all be absent.');
    end
    if ~all(cellfun(@isempty,axis_fields))
        axis_values=[axis_edge obs_freq proc_sweep];
        if (~isnumeric(axis_values))||(~isreal(axis_values))||...
           (numel(axis_values)~=3)||any(~isfinite(axis_values))||...
           (obs_freq<=0)||(proc_sweep<=0)
            error('OFFSET, SF, and SW_p must be real finite scalars with positive SF and SW_p.');
        end
    end
end
end

% All models are wrong, but some are useful.
%
% George E. P. Box


