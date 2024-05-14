% Saves pulses in Bruker format. The result is a text file with a 
% list of amplitudes and phases, usable in TopSpin. Syntax:
%
%                 bruker_write(X,Y,dt,file_name)
%
% Parameters:
%
%    X         - in-phase channel pulse amplitudes, a 
%                column vector in Hz
%
%    Y         - in-phase channel pulse amplitudes, a 
%                column vector in Hz
%
%    dt        - time slice duration, seconds
%
%    file_name - name of output file with .txt exten-
%                sion, a character string
%
% Outputs:
%
%    the function writes an ASCII text file
%
% a.acharya@soton.ac.uk
% david.goodwin@inano.au.dk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=bruker_write.m>

function bruker_write(X,Y,dt,file_name)

% Check consistency
grumble(X,Y,dt,file_name);

% Get amplitudes and phases
[amp,phi]=cartesian2polar(X,Y);

% Wrap phases and convert to degrees
phi=rad2deg(wrapTo2Pi(phi));

% Interval count, min and max
n_ints=numel(X); T=n_ints*dt;
min_amp=min(amp); max_amp=max(amp);
min_phi=min(phi); max_phi=max(phi);

% Date and time
date=datetime('now','Format','dd/MM/yyyy');
time=datetime('now','Format','HH:mm:ss');

% Write shaped pulse lines
lines = ["##TITLE= Optimal control pulse",  ... 
         "##JCAMP-DX= 5.00 Bruker JCAMP library",    ...
         "##DATA TYPE= Shape Data",                  ...
         "##ORIGIN= SPINACH",            ...
         "##OWNER= <demo>",                          ...
         "##DATE=" + string(date),                   ...
         "##TIME=" + string(time),                   ...
         "##$SHAPE_PARAMETERS=",                     ...
         "##MINX= " + num2str(min_amp),              ...
         "##MAXX= " + num2str(max_amp),              ...
         "##MINY= " + num2str(min_phi),              ...
         "##MAXY= " + num2str(max_phi),              ...
         "##$SHAPE_EXMODE= Universal",               ...
         "##$SHAPE_TOTROT=",                         ...
         "##$SHAPE_TYPE= Excitation",                ...
         "##$SHAPE_USER_DEF=",                       ...
         "##$SHAPE_REPHFAC= ",                       ...
         "##$SHAPE_BWFAC=",                          ...
         "##$SHAPE_BWFAC50=",                        ...
         "##$SHAPE_INTEGFAC=",                       ...
         "##$SHAPE_MODE=",                           ...
         "##$SHAPE_LENGTH=" + num2str(1e6*T),            ...
         "##NPOINTS="       + num2str(n_ints),       ...
         "##XYPOINTS= (XY..XY..)"];

% Append to output file 
writelines(lines,file_name);
writematrix([amp phi],file_name,'WriteMode','append','Delimiter',' ');
writelines('##END',file_name,'WriteMode','append');

end

% Consistency enforcement
function grumble(X,Y,dt,file_name)
if (~isnumeric(X))||(~isreal(X))||(~iscolumn(X))
    error('X must be a real column vector.');
end
if (~isnumeric(Y))||(~isreal(Y))||(~iscolumn(Y))
    error('X must be a real column vector.');
end
if numel(X)~=numel(Y)
    error('X and Y must have the same number of elements.');
end
if (~isnumeric(dt))||(~isreal(dt))||...
   (~isscalar(dt))||(dt<=0)
    error('dt must be a positive real scalar.');
end
if ~ischar(file_name)
    error('file_name must be a character string.');
end
end

% We are stuck with technology when all we really 
% want is just stuff that works. How do you recog-
% nise something that is still technology? A good
% clue is if it comes with a manual.
%
% Douglas Adams

