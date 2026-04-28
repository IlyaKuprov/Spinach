% Spin system of 13C-labelled lactate with the OH protons
% assumed to be in rapid exchange with water. Syntax:
%
%             [sys,inter,bas]=lactate(spins)
%
% Parameters:
%
%    spins - which spins to include, e.g. {'1H','13C'}
%
% Outputs:
%
%    sys   - Spinach spin system description structure
%
%    inter - Spinach interaction description structure
%
% elton.montrazi@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=lactate.m>

function [sys,inter]=lactate(spins)

% Check consistency
grumble(spins);

% Isotopes
sys.isotopes={'13C','13C','13C','1H','1H','1H','1H'};
sys.labels={'CO','CA','CB','HA','HB1','HB2','HB3'};

% Chemical shifts (approximate)
inter.zeeman.scalar={182.0 68.4 20.6 4.0 1.2 1.2 1.2};

% J-couplings (very approximate)
inter.coupling.scalar=cell(7,7);
inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'HA')}=3.9;
inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'HB1')}=3.7;
inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'HB2')}=3.7;
inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'HB3')}=3.7;
inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'CA')}=40.0;
inter.coupling.scalar{idxof(sys,'CO'),idxof(sys,'CB')}=3.0;
inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'HA')}=145.0;
inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'HB1')}=4.4;
inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'HB2')}=4.4;
inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'HB3')}=4.4;
inter.coupling.scalar{idxof(sys,'CA'),idxof(sys,'CB')}=40.0;
inter.coupling.scalar{idxof(sys,'CB'),idxof(sys,'HB1')}=128.0;
inter.coupling.scalar{idxof(sys,'CB'),idxof(sys,'HB2')}=128.0;
inter.coupling.scalar{idxof(sys,'CB'),idxof(sys,'HB3')}=128.0;
inter.coupling.scalar{idxof(sys,'CB'),idxof(sys,'HA')}=3.6;
inter.coupling.scalar{idxof(sys,'HA'),idxof(sys,'HB1')}=6.9;
inter.coupling.scalar{idxof(sys,'HA'),idxof(sys,'HB2')}=6.9;
inter.coupling.scalar{idxof(sys,'HA'),idxof(sys,'HB3')}=6.9;

% Prune the arrays 
mask=ismember(sys.isotopes,spins);
sys.isotopes=sys.isotopes(mask);
sys.labels=sys.labels(mask);
inter.zeeman.scalar=inter.zeeman.scalar(mask);
inter.coupling.scalar=inter.coupling.scalar(mask,mask);

end

% Consistency enforcement
function grumble(spins)
if (~iscell(spins))||(~all(cellfun(@ischar,spins)))
    error('spins argument must be a cell array of character strings.');
end
end

% Why 100? If I were wrong, one would have been enough.
%
% Albert Eistein's response to Nazis lining up
% 100 Aryan scientists to denounce his theories

