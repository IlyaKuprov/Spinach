% Relaxation rate of Gd(III) as a function of zero-field
% splitting, computed using Redfield theory. The correla-
% tion time (1 fs) refers to vibrational dynamics of the
% ligand cage.
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk
% e.suturina@soton.ac.uk

function lanthanide_redfield()

% Spin system properties
sys.isotopes={'E8'};
inter.zeeman.scalar={1.9918};

% Magnet field
sys.magnet=9.40;

% Relaxation parameters
inter.relaxation={'redfield'};
inter.rlx_keep='labframe';
inter.equilibrium='zero';
inter.tau_c={1e-15};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Linearly spaced B20 in cm^-1
b20=linspace(0.1,10,10)';

% Loop over B20 values
for n=1:numel(b20)
    
    % Giant spin Hamiltonian parameters
    inter.giant.coeff={{[0 0 0],[0 0 icm2hz(b20(n)) 0 0]}};
    inter.giant.euler={{[0 0 0],[0 0 0]}};
    
    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Relaxation superoperator
    R=relaxation(spin_system);

    % Longitudinal and tranverse relaxation rates
    Lp=state(spin_system,'L+','E8'); 
    Lz=state(spin_system,'Lz','E8');
    T1(n)=-1/((Lz'*R*Lz)/(Lz'*Lz)); %#ok<AGROW>
    T2(n)=-1/((Lp'*R*Lp)/(Lp'*Lp)); %#ok<AGROW>

end

% Do the plotting
figure(); plot(b20,[T1; T2],'o');
set(gca,'yscale','log'); kgrid;
kxlabel('ZFS $B_{2}^{0} (cm^{-1})$');
kylabel('Relaxation time, seconds');

end

