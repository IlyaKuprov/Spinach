% Resets the stupid ass figure defaults in R2025a 
% and later back to sensible values.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=pauli.m>

function kfigure(varargin)

% Reset to pre-R2025a settings
set(groot,'defaultFigurePosition',[680 458 560 420]); 
set(groot,'defaultFigureWindowStyle','normal'); 
set(groot,'defaultFigureMenuBar','figure'); 
set(groot,'defaultFigureToolbar','figure'); 

% Create
figure(varargin{:});

end

% #NGRUM #NHEAD

% The most common error of a smart engineer is to 
% optimize a thing that should not exist.
%
% Elon Musk

