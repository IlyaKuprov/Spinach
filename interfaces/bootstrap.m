% A minimal spin_system structure required to call many
% Spinach functions. Use this function if your system is 
% not a spin system, and you simply want to use one of
% Spinach functions outside the context. Syntax:
%
%             spin_system=bootstrap(volume)
%
% Parameters:
%
%    volume - the content of sys.output in the bootstrap
%             call, set to 'hush' to make the resulting
%             object suppress console output
%
% Outpus:
%
%    spin_system  - empty but valid Spinach spin system
%                   description object
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bootstrap.m>

function spin_system=bootstrap(volume)

% Default volume
if ~exist('volume','var')
    volume='console';
end

% Check consistency
grumble(volume);

% Ghost spin
sys.magnet=0; sys.isotopes={'G'};

% No spin interactions
inter.zeeman.matrix=cell(1);
inter.coupling.matrix=cell(1);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
if exist('volume','var')
    sys.output=volume;
end
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

end

% Consistency enforcement
function grumble(volume)
if ~ischar(volume)
    error('volume must be a character string.');
end
end

% "Without love, humans would be... rare."
%
% Terry Pratchett

