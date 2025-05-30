% Transverse relaxation rate as a function of the applied magnetic 
% field in a 3-fluorotyrosine labelled protein. The fluorine atom
% and its directly bonded carbon are included.
% 
% Calculation time: minutes.
%
% ilya.kuprov@weizmann.ac.il

function trosy_fluorine_num()

% Read 3-fluorotyrosine DFT calculation
[~,inter_dft]=g2spinach(gparse('../standard_systems/3_fluoro_tyr.log'),...
                            {{'C','13C'},{'F','19F'}},[186.38 192.97]);
                                 
% Extract coordinates and CSAs
sys.isotopes={'19F','13C'};
inter.zeeman.matrix=cell(1,2);
inter.zeeman.matrix{1}=inter_dft.zeeman.matrix{8};
inter.zeeman.matrix{2}=inter_dft.zeeman.matrix{7};
inter.coordinates=cell(2,1);
inter.coordinates{1}=inter_dft.coordinates{8};
inter.coordinates{2}=inter_dft.coordinates{7};

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
    LpF=state(spin_system,{'L+'},{1});
    LpC=state(spin_system,{'L+'},{2});
    LpF=LpF/norm(LpF,2); LpC=LpC/norm(LpC,2);
    F_left=state(spin_system,{'L+'},{1})-...
         2*state(spin_system,{'L+','Lz'},{1,2});
    F_right=state(spin_system,{'L+'},{1})+...
          2*state(spin_system,{'L+','Lz'},{1,2});
    C_left=state(spin_system,{'L+'},{2})-...
         2*state(spin_system,{'L+','Lz'},{2,1});
    C_right=state(spin_system,{'L+'},{2})+...
          2*state(spin_system,{'L+','Lz'},{2,1});
    F_left=F_left/norm(F_left,2);
    F_right=F_right/norm(F_right,2);
    C_left=C_left/norm(C_left,2);
    C_right=C_right/norm(C_right,2);
    
    % Relaxation rates
    r2c(n)=-LpC'*R*LpC; 
    r2f(n)=-LpF'*R*LpF;
    fleft(n)=-F_left'*R*F_left;
    fright(n)=-F_right'*R*F_right;
    cleft(n)=-C_left'*R*C_left;
    cright(n)=-C_right'*R*C_right;
    
end
                         
% Plotting
figure();
plot(lin_freq',[fleft' r2f' fright']); 
kxlabel('Proton Larmor frequency, MHz'); kgrid;
kylabel('Relaxation matrix element, Hz'); xlim tight;
klegend({'${\hat F_ + } - 2{\hat F_ + }{\hat C_{\rm{Z}}}$',...
         '${\hat F_ + }$',...
         '${\hat F_ + } + 2{\hat F_ + }{\hat C_{\rm{Z}}}$'},...
         'Location','northwest','FontSize',12);
    
figure();
plot(lin_freq',[cleft' r2c' cright']); 
kxlabel('Proton Larmor frequency, MHz'); kgrid;
kylabel('Relaxation matrix element, Hz'); xlim tight;
klegend({'${\hat C_ + } - 2{\hat C_ + }{\hat F_{\rm{Z}}}$',...
         '${\hat C_ + }$',...
         '${\hat C_ + } + 2{\hat C_ + }{\hat F_{\rm{Z}}}$'},...
         'Location','northwest','FontSize',12);

end

                         