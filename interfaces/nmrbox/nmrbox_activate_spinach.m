% Spinach activation script for NMRBox.
%
% ilya.kuprov@weizmann.ac.uk
%
% #NGRUM #NWIKI #NHEAD

% Set up paths
addpath(genpath('/usr/software/spinach/etc'));
addpath(genpath('/usr/software/spinach/experiments'));
addpath(genpath('/usr/software/spinach/interfaces'));
addpath(genpath('/usr/software/spinach/kernel'));

% Test and display the welcome message
existentials(); disp('Welcome to Spinach 2.9');
disp('Type ''wiki'' and press Enter to open documentation web site.');
disp('Type ''forum'' and press Enter to open user support forum.');
disp('See /usr/software/spinach/examples directory for the example set.');

% Go into the examples directory
cd('/usr/software/spinach/examples');

% Sir Isaac Newton was rigorously puritanical: when one of
% his few friends told him "a loose story about a nun", he
% ended their friendship. He is not known to have ever had
% a romantic relationship of any kind, and is believed to
% have died a virgin.

