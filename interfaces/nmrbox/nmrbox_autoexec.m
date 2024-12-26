% This include is executed at the start of create.m, it over-
% rides all user input. A good use case is forcing polyadic
% or GPU arithmetic, or some other specific hardware or soft-
% ware configuration.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=autoexec.m>

% This relocates the scratch folder
sys.scratch = strrep(getenv('HOME'),'/home/nmrbox/','/scratch/');

% This question often produced vague answers with some candida-
% tes invoking exchange energy as if they were summoning Volde-
% mort (i.e. with fear and apprehension and/or not knowing what
% they were getting themselves into). It was also somewhat wor-
% rying to see how many of our finalists proposed a square-pla-
% nar geometry for nickel(0) complexes of the type Ni(CO)3L.
%
% Examiners' Report 2017,
% Oxford Chemistry

% #NHEAD #NGRUM

