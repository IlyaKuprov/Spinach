% Transverse relaxation rate as a function of the applied magnetic 
% field at a typical amide N-H group in a protein. Rotational cor-
% relation time set to 25 ns. Nitrogen CSA parameters from
%
%              http://dx.doi.org/10.1021/ja0016194
%
% Nitrogen-proton bond length from DFT.
%
% Calculation time: minutes.
% 
% ilya.kuprov@weizmann.ac.il

function trosy_nh()

% Specify coordinates and CSAs
sys.isotopes={'1H','15N'};
inter.zeeman.eigs{1}=[   6   0   -6];
inter.zeeman.eigs{2}=[-108   62  46];
inter.zeeman.euler{1}=pi*[0   0   0]/180;
inter.zeeman.euler{2}=pi*[0   0 -19]/180;
inter.coordinates{1}=[1.04 0.00 0.00];
inter.coordinates{2}=[0.00 0.00 0.00];

% Relaxation theory
inter.relaxation={'redfield'};
inter.rlx_keep='labframe';
inter.equilibrium='zero';
inter.tau_c={25e-9};

% Formalism and approximation
bas.formalism='sphten-liouv';
bas.approximation='none';

% Disable startup checks
sys.disable={'hygiene'};

% Magnetic field grid
lin_freq=linspace(200,1500,30);
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
    LpN=state(spin_system,{'L+'},{2});
    LpH=LpH/norm(LpH,2); LpN=LpN/norm(LpN,2);
    H_left=state(spin_system,{'L+'},{1})-...
         2*state(spin_system,{'L+','Lz'},{1,2});
    H_right=state(spin_system,{'L+'},{1})+...
          2*state(spin_system,{'L+','Lz'},{1,2});
    N_left=state(spin_system,{'L+'},{2})-...
         2*state(spin_system,{'L+','Lz'},{2,1});
    N_right=state(spin_system,{'L+'},{2})+...
          2*state(spin_system,{'L+','Lz'},{2,1});
    H_left=H_left/norm(H_left,2);
    H_right=H_right/norm(H_right,2);
    N_left=N_left/norm(N_left,2);
    N_right=N_right/norm(N_right,2);
    
    % Relaxation rates
    r2c(n)=-LpN'*R*LpN; 
    r2h(n)=-LpH'*R*LpH;
    hleft(n)=-H_left'*R*H_left;
    hright(n)=-H_right'*R*H_right;
    nleft(n)=-N_left'*R*N_left;
    nright(n)=-N_right'*R*N_right;
    
end
                         
% Plotting
figure(); plot(lin_freq',[hleft' r2h' hright']);
kxlabel('Proton Larmor frequency, MHz'); kgrid;
kylabel('Relaxation matrix element, Hz'); xlim tight;
klegend({'${\hat H_ + } - 2{\hat H_ + }{\hat N_{\rm{Z}}}$',...
         '${\hat H_ + }$',...
         '${\hat H_ + } + 2{\hat H_ + }{\hat N_{\rm{Z}}}$'},...
         'Location','northwest','FontSize',12);
    
figure(); plot(lin_freq',[nleft' r2c' nright']);
kxlabel('Proton Larmor frequency, MHz'); kgrid;
kylabel('Relaxation matrix element, Hz'); xlim tight;
klegend({'${\hat N_ + } - 2{\hat N_ + }{\hat H_{\rm{Z}}}$',...
         '${\hat N_ + }$',...
         '${\hat N_ + } + 2{\hat N_ + }{\hat H_{\rm{Z}}}$'},...
         'Location','northwest','FontSize',12);

end

                         