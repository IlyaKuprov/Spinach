% A replacement for the 'grid' command in Matlab that
% produces grey (rather than black-and-transparent) grid
% lines that are suitable for publishing.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=kgrid.m>

function kgrid()

% Publisher-friendly grid settings
grid on; grid minor; 
set(gca,'MinorGridAlpha',1,...
        'MinorGridColor',[0.95 0.95 0.95],...
        'MinorGridLineStyle','-')
set(gca,'GridAlpha',1,...
        'GridColor',[0.85 0.85 0.85],...
        'GridLineStyle','-');
set(gca,'Layer','bottom');

end

% Frigid gentlewomen of the jury! I had thought that 
% months, perhaps years, would elapse before I dared
% to reveal myself to Dolores Haze; but by six she 
% was wide awake, and by six fifteen we were techni-
% cally lovers. I am going to tell you something very
% strange: it was she who seduced me.
%
% Vladimir Nabokov, "Lolita"

% #NGRUM #NHEAD