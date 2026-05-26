% Spin system generator for directional derivative tests. Syntax:
%
%    [spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalism)
%
% ilya.kuprov@weizmann.ac.il

function [spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalism)

% Check consistency
grumble(formalism);

% Select system size
switch formalism

    case 'sphten-liouv'

        % Keep the original large Liouville-space test
        n_spins=100;

    case {'zeeman-liouv','zeeman-hilb'}

        % Use a compact system for full Zeeman formalisms
        n_spins=2;

end

% Set the magnetic field
sys.magnet=28.18;
sys.output='hush';

% Put non-interacting spins at equal intervals
% within the [-100,+100] ppm chemical shift range
sys.isotopes=cell(n_spins,1);
for n=1:n_spins
    sys.isotopes{n}='13C';
end
inter.zeeman.scalar=num2cell(linspace(-100,100,n_spins));

% Select the requested basis set
bas.formalism=formalism;
switch formalism

    case 'sphten-liouv'

        % Keep complete single-spin terms only
        bas.approximation='IK-2';
        bas.space_level=1;
        bas.connectivity='scalar_couplings';

    case {'zeeman-liouv','zeeman-hilb'}

        % Keep the full Zeeman basis
        bas.approximation='none';

end

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up spin states
Sx=state(spin_system,'Lx','13C'); Sx=Sx/norm(full(Sx),2);
Sy=state(spin_system,'Ly','13C'); Sy=Sy/norm(full(Sy),2);
Sz=state(spin_system,'Lz','13C'); Sz=Sz/norm(full(Sz),2);

% Get control operators
Lx=operator(spin_system,'Lx','13C');
Ly=operator(spin_system,'Ly','13C');

% Get the drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

end

% Consistency enforcement
function grumble(formalism)
if (~ischar(formalism))||(~ismember(formalism,{'sphten-liouv','zeeman-liouv','zeeman-hilb'}))
    error('formalism must be ''sphten-liouv'', ''zeeman-liouv'', or ''zeeman-hilb''.');
end
end
