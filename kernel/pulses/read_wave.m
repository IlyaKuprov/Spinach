% Reads JCAMP-DX pulse waveform files (a few examples are distri-
% buted with Spinach, see /kernel/pulses/pk_files). Syntax:
%
%   [A,phi,Cx,Cy,scaling_factor]=read_wave(filename,npoints)
%
% Parameters:
%
%   filename   -   a string containing the name of the file
%
%   npoints    -   waveform upsampling or downsampling is
%                  performed to this number of points
%
% Outputs:
%
%   A              -   polar amplitude at each slice
%
%   phi            -   polar phase at each slice, radians
%
%   Cx             -   Cartesian amplitude in X at each slice
%
%   Cy             -   Cartesian amplitude in Y at each slice
%
%   scaling factor -   scaling factor for a given pulse shape
% 
% Note: put your own pulses into /kernel/pulses/pk_files; please
%       also consider sending them to us.
%
% ilya.kuprov@weizmann.ac.il
% kpervushin@ntu.edu.sg
% mariagrazia.concilio@sjtu.edu.cn
%
% <https://spindynamics.org/wiki/index.php?title=read_wave.m>

function [A,phi,Cx,Cy,scaling_factor]=read_wave(filename,npoints)

% Check consistency
grumble(filename,npoints);

% Read the waveform
P=mfilename('fullpath'); P=P(1:(end-9));
wavefile=fopen([P 'pk_files/' filename],'r');
waveform=textscan(wavefile,'%f, %f','CommentStyle','##'); 
waveform=cell2mat(waveform); frewind(wavefile);

% Read the scaling factor
line_by_line=textscan(wavefile,'%s','delimiter','\n');
fclose(wavefile); line_by_line=line_by_line{1};
for n=1:numel(line_by_line)
    if (numel(line_by_line{n})>17)&&...
        strcmp(line_by_line{n}(1:17),'##$SHAPE_INTEGFAC')
        scaling_factor=str2double(line_by_line{n}(19:end));
    end
end
if (~exist('scaling_factor','var'))||isnan(scaling_factor)
    error('the scaling factor string '' ##$SHAPE_INTEGFAC= '' must be present in the pk_file.');
end

% Get amplitude and phase
A=waveform(:,1)/100; phi=pi*waveform(:,2)/180;

% Resample amplitude and phase
A=interp1(linspace(0,1,numel(A)),A,linspace(0,1,npoints),'pchip');
phi=interp1(linspace(0,1,numel(phi)),phi,linspace(0,1,npoints),'pchip');

% Convert into Cartesian coordinates
if nargout>3,[Cx,Cy]=polar2cartesian(A,phi); end

end

% Consistency enforcement
function grumble(filename,npoints)
if ~ischar(filename)
    error('filename must be a character string.');
end
if (~isnumeric(npoints))||(~isreal(npoints))||(~isfinite(npoints))||...
   (numel(npoints)~=1)||(npoints<1)||mod(npoints,1)
    error('npoints must be a positive real integer.');
end
end

% I condemn Christianity... It is, to me, the greatest of all imaginable
% corruptions. [...] To breed in humans a self-contradiction, an art of
% self-pollution, a will to lie at any price, an aversion and contempt for
% all good and honest instincts! [...] The beyond as the will to deny all
% reality; the cross as the distinguishing mark of the most subterranean
% conspiracy ever heard of -- against health, beauty, well-being, intel-
% lect... against life itself.
%
% I call Christianity the one great curse, the one great intrinsic depra-
% vity, the one great instinct of revenge, for which no means are venomo-
% us enough, or secret, subterranean and small enough - I call it the one
% immortal blemish upon the human race...
%
% Friedrich Nietzsche

