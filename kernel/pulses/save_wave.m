% Writes a pulse waveform to a JCAMP-DX pulse file. Syntax:
%
%              save_wave(filename,Cx,Cy,npoints,param)
%
% Parameters:
%
%    filename  -   a string containing the name of the file
%
%    Cx        -   Cartesian amplitude in X at each slice, normalised
%
%    Cy        -   Cartesian amplitude in Y at each slice, normalised
%
%    npoints   -   waveform upsampling or downsampling is performed
%                  to this number of points.
%
%    param.acc         - [optional] format string controlling
%                        output accuracy, default is '%.6e'
%
%    param.pwr_level   - pulse power, rad/s
% 
%    param.pulse_dt    - a vector of pulse slice durations, seconds
% 
%    param.bandwidth   - pulse bandwidth, Hz
%
%    param.alpha       - pulse flip angle, radians.
%
%    param.title       - [optional] JCAMP title string
%
%    param.owner       - [optional] JCAMP owner string
%
%    param.shape_mode  - JCAMP shape mode string
%
%    param.shape_type  - JCAMP shape type string
%
%    param.shape_param - [optional] JCAMP shape parameter string
% 
% Outputs:
%
%    the function writes an ASCII text file
%
% david.goodwin@chem.ox.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=save_wave.m>

function save_wave(filename,Cx,Cy,npoints,param)

% Check consistency
grumble(filename,Cx,Cy,npoints,param);

% Set default accuracy
if ~isfield(param,'acc'), param.acc='%.6e'; end

% Resample waveform
Cx=interp1(linspace(0,1,numel(Cx)),Cx,linspace(0,1,npoints),'pchip');
Cy=interp1(linspace(0,1,numel(Cy)),Cy,linspace(0,1,npoints),'pchip');

% Convert into phase-amplitude representation
[Ca,Cp]=cartesian2polar(param.pwr_level*Cx,param.pwr_level*Cy);

% Stretch, normalise, and wrap phases
Ca=100*Ca(:)/param.pwr_level; Cp=wrapTo360(rad2deg((Cp(:))));

% Clean up the numerics
Ca=1e-12*round((1/1e-12)*Ca); Cp=1e-12*round((1/1e-12)*Cp);

% Make sure amplitudes are reasonable
if any(Ca>100), error('waveform amplitude exceeds param.pwr_level'); end

% Calculate pulse duration
pulse_dur=sum(param.pulse_dt);

% Find maximum and minimum values
min_amp=pk_strfmt(min(Ca),param.acc);
max_amp=pk_strfmt(max(Ca),param.acc);
min_phi=pk_strfmt(min(Cp),param.acc);
max_phi=pk_strfmt(max(Cp),param.acc);

% Date and time
date_string=datestr(clock,'yyyy/mm/dd'); %#ok<CLOCK,DATST> 
time_string=datestr(clock,'HH:MM:SS');   %#ok<CLOCK,DATST> 

% Bandwidth and scaling factors
bandwidth_factor=pulse_dur*param.bandwidth;
bandwidth_factor=pk_strfmt(bandwidth_factor,param.acc);
scaling_factor=param.alpha/(param.pwr_level*pulse_dur);
scaling_factor=pk_strfmt(scaling_factor,param.acc);

% Rotation angle
alpha=pk_strfmt(wrapTo360(rad2deg(param.alpha)),param.acc);

% Defaults for the header
if ~isfield(param,'title'), param.title=''; end
if ~isfield(param,'owner'), param.owner=''; end
if ~isfield(param,'shape_param'), param.shape_param=''; end
if ~isfield(param,'spinach_ver'), param.spinach_ver='SPINACH'; end

% Open the file for writing
wavefile=fopen(filename,'w+');

% Write the header
fprintf(wavefile,...
        ['##TITLE= ' param.title '\n',...
         '##JCAMP-DX= 5.00 Bruker JCAMP library\n',...
         '##DATA TYPE= Shape Data\n',...
         '##ORIGIN= ' param.spinach_ver  '\n',...
         '##OWNER= ' param.owner '\n',...
         '##DATE= ' date_string '\n',...
         '##TIME= ' time_string '\n',...
         '##$SHAPE_PARAMETERS= Type: ' param.shape_param '\n',...
         '##MINX= ' min_amp '\n',...
         '##MAXX= ' max_amp '\n',...
         '##MINY= ' min_phi '\n',...
         '##MAXY= ' max_phi '\n',...
         '##$SHAPE_EXMODE= ' param.shape_mode '\n',...
         '##$SHAPE_TOTROT= ' alpha '\n',...
         '##$SHAPE_TYPE= ' param.shape_type '\n',...
         '##$SHAPE_USER_DEF= \n',...
         '##$SHAPE_REPHFAC= \n',...
         '##$SHAPE_BWFAC= ' bandwidth_factor '\n',...
         '##$SHAPE_BWFAC50= \n',...
         '##$SHAPE_INTEGFAC= ' scaling_factor '\n',...
         '##$SHAPE_MODE= 0\n',...
         '##NPOINTS= ' num2str(npoints,'%.0f') '\n',...
         '##XYPOINTS= (XY..XY)\n']);
    
% Write the waveform
waveform_string=sprintf([param.acc ', ' param.acc '\n'],[Ca(:) Cp(:)]');
fprintf(wavefile,pk_strfmt(waveform_string,param.acc));

% Write the terminator and close the file
fprintf(wavefile,'##END'); fclose(wavefile);

end

% Convert numbers to strings in JCAMP-DX format
function str_out=pk_strfmt(num_in,acc)

% Format in standard form
str_out=num2str(num_in,acc);

% Replace lowercase e with uppercase
str_out(str_out=='e')='E';

% Remove + in exponent
str_out(str_out=='+')='';

end

% Consistency enforcement
function grumble(filename,Cx,Cy,npoints,param)
if ~ischar(filename)
    error('filename must be a character string.');
end
if (~isnumeric(Cx))||(~isreal(Cx))||(~isrow(Cx))
    error('Cx must be a row vector of real numbers.');
end
if (~isnumeric(Cy))||(~isreal(Cy))||(~isrow(Cy))
    error('Cy must be a row vector of real numbers.');
end
if length(Cx)~=length(Cy)
    error('Cx and Cy must be vectors of the same length.')
end
if ~isfield(param,'pulse_dt')
    error('param.pulse_dt is required')
elseif (~isnumeric(param.pulse_dt))||(~isreal(param.pulse_dt))||...
     any(~isfinite(param.pulse_dt))||any(param.pulse_dt<0)
    error('param.pulse_dt must be a vector of non-negative real numbers.');
elseif length(Cx)~=length(param.pulse_dt)
    error('param.pulse_dt must be a vector of the same length as Cx and Cy.')
end
if (~isnumeric(npoints))||(~isreal(npoints))||(~isfinite(npoints))||...
   (numel(npoints)~=1)||(npoints<1)||mod(npoints,1)
    error('npoints must be a positive real integer.');
end
if ~isfield(param,'bandwidth')
    error('param.bandwidth is required')
elseif (~isnumeric(param.bandwidth))||(~isreal(param.bandwidth))||...
        (numel(param.bandwidth)~=1)||(~isfinite(param.bandwidth))
    error('param.bandwidth must be a real number.')
end
if ~isfield(param,'pwr_level')
    error('param.pwr_level is required')
elseif (~isnumeric(param.pwr_level))||(~isreal(param.pwr_level))||...
        (numel(param.pwr_level)~=1)||(~isfinite(param.pwr_level))
    error('param.pwr_level must be a real number.')
end
if ~isfield(param,'alpha')
    error('param.alpha is required')
elseif (~isnumeric(param.alpha))||(~isreal(param.alpha))||...
        (numel(param.alpha)~=1)||(~isfinite(param.alpha))
    error('param.alpha must be a real number.')
end
if isfield(param,'title')&&~ischar(param.title)
    error('param.title must be a character string.')
end
if isfield(param,'owner')&&~ischar(param.owner)
    error('param.owner must be a character string.')
end
if isfield(param,'shape_mode')&&~ischar(param.shape_mode)
    error('param.shape_mode must be a character string.')
end
if isfield(param,'shape_type')&&~ischar(param.shape_type)
    error('param.shape_type must be a character string.')
end
if isfield(param,'shape_param')&&~ischar(param.shape_param)
    error('param.shape_param must be a character string.')
end
end

% Nothing was your own except the few
% cubic centimetres inside your skull.
%
% George Orwell, "1984"

