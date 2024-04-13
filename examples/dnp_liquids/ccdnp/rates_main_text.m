% Self- and cross-relaxation rates in cross-correlated DNP, considering
% a system with two electrons connected by exchange coupling, both cou-
% pled to a nucleus by dipolar couplings. Further particulars in:
%
%               https://doi.org/10.1016/j.jmr.2021.106940
%
% Calculation time: seconds
%
% m.g.concilio@soton.ac.uk 
% i.kuprov@soton.ac.uk

function rates_main_text()

% Magnet field
sys.magnet=14.1;

% Spin system
sys.isotopes={'1H','E','E'};
  
% Zeeman interactions
inter.zeeman.eigs{1,1}=[0 10 20];
inter.zeeman.euler{1,1}=[0 0 0];   
inter.zeeman.eigs{1,2}=[2.003400 2.003800 2.003800];          
inter.zeeman.euler{1,2}=[-0.872 -0.013 0.868];     
inter.zeeman.eigs{1,3}=[2.005700 2.003000 2.003000];  
inter.zeeman.euler{1,3}=[-1.145 0.061 1.143];

% Exchange coupling
inter.coupling.scalar=cell(3,3);
inter.coupling.scalar{2,3}=3e6;        

% Cooridnates
inter.coordinates={[ 0.000  0.000  0.000];
                   [ 5.090  0.010  0.958];
                   [-5.090  0.061  1.032]};   
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.temperature=298;
inter.tau_c={100e-12};          
sys.tols.rlx_integration=1e-10;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relaxation superoperator
R=relaxation(spin_system);

% Rate printouts
disp(' ');
Nz=state(spin_system,{'Lz'},{1});  Nz=Nz/norm(Nz,2);                
E1z=state(spin_system,{'Lz'},{2}); E1z=E1z/norm(E1z,2);                
E2z=state(spin_system,{'Lz'},{3}); E2z=E2z/norm(E2z,2);
disp(['R(E1z): ' num2str(E1z'*R*E1z)]);
disp(['R(E2z): ' num2str(E2z'*R*E2z)]);
disp(['R(E1z -> Nz): ' num2str(E1z'*R*Nz)]);
disp(['R(E2z -> Nz): ' num2str(E2z'*R*Nz)]);
disp(['R(Nz):  ' num2str(Nz'*R*Nz)]);

disp(' ');
E1p=state(spin_system,{'L+'},{2});              
E2p=state(spin_system,{'L+'},{3});
E1pE2z=state(spin_system,{'L+','Lz'},{2,3});
E1zE2p=state(spin_system,{'Lz','L+'},{2,3});
E1_a=E1p+2*E1pE2z; E1_a=E1_a/norm(E1_a,2);
E1_b=E1p-2*E1pE2z; E1_b=E1_b/norm(E1_b,2);
E2_a=E2p+2*E1zE2p; E2_a=E2_a/norm(E2_a,2);
E2_b=E2p-2*E1zE2p; E2_b=E2_b/norm(E2_b,2);
disp(['R(E1p+2*E1pE2z): ' num2str(E1_a'*R*E1_a)]);
disp(['R(E1p-2*E1pE2z): ' num2str(E1_b'*R*E1_b)]);
disp(['R(E2p+2*E1zE2p): ' num2str(E2_a'*R*E2_a)]);
disp(['R(E2p-2*E1zE2p): ' num2str(E2_b'*R*E2_b)]);

disp(' ');
E1p=state(spin_system,{'L+'},{2}); E1p=E1p/norm(E1p,2);      
NzE1p=state(spin_system,{'Lz','L+'},{1,2}); NzE1p=NzE1p/norm(NzE1p,2);
E1pE2z=state(spin_system,{'L+','Lz'},{2,3}); E1pE2z=E1pE2z/norm(E1pE2z,2);
NzE1pE2z=state(spin_system,{'Lz','L+','Lz'},{1,2,3}); NzE1pE2z=NzE1pE2z/norm(NzE1pE2z,2);
disp(['R([E1p -> NzE1p]): ' num2str(NzE1p'*R*E1p)]);
disp(['R([E1p -> NzE1pE2z]): ' num2str(NzE1pE2z'*R*E1p)]);
disp(['R([E1pE2z -> NzE1pE2z]): ' num2str(NzE1pE2z'*R*E1pE2z)]);

disp(' ');
Nz=state(spin_system,{'Lz'},{1}); Nz=Nz/norm(Nz,2); 
NzE1z=state(spin_system,{'Lz','Lz'},{1,2}); NzE1z=NzE1z/norm(NzE1z,2);         
NzE2z=state(spin_system,{'Lz','Lz'},{1,3}); NzE2z=NzE2z/norm(NzE2z,2);
NzE1zE2z=state(spin_system,{'Lz','Lz','Lz'},{1,2,3}); NzE1zE2z=NzE1zE2z/norm(NzE1zE2z,2);
disp(['R(NzE1z -> Nz): ' num2str(Nz'*R*NzE1z)]);
disp(['R(NzE2z -> Nz): ' num2str(Nz'*R*NzE2z)]);
disp(['R(NzE1zE2z -> Nz): ' num2str(Nz'*R*NzE1zE2z)]);

end

