% Transverse relaxation rate as a function of the applied magnetic field
% at the C-H group in position 3 of the aromatic ring of tyrosine.
%
% Calculation time: minutes.
% 
% ilya.kuprov@weizmann.ac.il

function trosy_proton()

% Read 3-fluorotyrosine DFT calculation
[~,inter_dft]=g2spinach(gparse('../standard_systems/amino_acids/tyr.log'),...
                                 {{'C','13C'},{'H','1H'}},[186.38 33.44]);
                                 
% Extract coordinates and CSAs
sys.isotopes={'1H','13C'};
inter.zeeman.matrix=cell(1,2);
inter.zeeman.matrix{1}=inter_dft.zeeman.matrix{10};
inter.zeeman.matrix{2}=inter_dft.zeeman.matrix{5};
inter.coordinates=cell(2,1);
inter.coordinates{1}=inter_dft.coordinates{10};
inter.coordinates{2}=inter_dft.coordinates{5};

% Relaxation theory
inter.relaxation={'redfield'};
inter.rlx_keep='labframe';
inter.equilibrium='zero';
inter.tau_c={25e-9};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Disable startup checks
sys.disable={'hygiene'};

% Magnetic field grid
lin_freq=linspace(200,800,20);
B0=2*pi*lin_freq*1e6/spin('1H');

% Loop over magnetic fields
for n=1:numel(B0) %#ok<*AGROW>
    
    % Set the magnet field
    sys.magnet=B0(n);
    
    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    
    % Relaxation superoperator
    R=relaxation(spin_system);
    
    % States of interest
    LpH=state(spin_system,{'L+'},{1});
    LpC=state(spin_system,{'L+'},{2});
    LpH=LpH/norm(LpH,2); LpC=LpC/norm(LpC,2);
    H_left=state(spin_system,{'L+'},{1})-...
         2*state(spin_system,{'L+','Lz'},{1,2});
    H_right=state(spin_system,{'L+'},{1})+...
          2*state(spin_system,{'L+','Lz'},{1,2});
    C_left=state(spin_system,{'L+'},{2})-...
         2*state(spin_system,{'L+','Lz'},{2,1});
    C_right=state(spin_system,{'L+'},{2})+...
          2*state(spin_system,{'L+','Lz'},{2,1});
    H_left=H_left/norm(H_left,2);
    H_right=H_right/norm(H_right,2);
    C_left=C_left/norm(C_left,2);
    C_right=C_right/norm(C_right,2);
    
    % Relaxation rates
    r2c(n)=-LpC'*R*LpC; 
    r2h(n)=-LpH'*R*LpH;
    hleft(n)=-H_left'*R*H_left;
    hright(n)=-H_right'*R*H_right;
    cleft(n)=-C_left'*R*C_left;
    cright(n)=-C_right'*R*C_right;
    
end
                         
% Plotting
figure();
plot(lin_freq',[hleft' r2h' hright']); ylim([0 160]);
kxlabel('Proton Larmor frequency, MHz'); kgrid;
kylabel('Relaxation matrix element, Hz');
klegend({'${\hat H_ + } - 2{\hat H_ + }{\hat C_{\rm{Z}}}$',...
         '${\hat H_ + }$',...
         '${\hat H_ + } + 2{\hat H_ + }{\hat C_{\rm{Z}}}$'},...
         'Location','northwest','FontSize',12);
    
figure();
plot(lin_freq',[cleft' r2c' cright']); ylim([0 400]);
kxlabel('Proton Larmor frequency, MHz'); kgrid;
kylabel('Relaxation matrix element, Hz');
klegend({'${\hat C_ + } - 2{\hat C_ + }{\hat H_{\rm{Z}}}$',...
         '${\hat C_ + }$',...
         '${\hat C_ + } + 2{\hat C_ + }{\hat H_{\rm{Z}}}$'},...
         'Location','northwest','FontSize',12);

end

                         