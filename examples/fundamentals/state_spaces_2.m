% Contributions from different orders of spin correlation to the system 
% trajectory in the pulse-acquire 1H NMR simulation of anti-3,5-difluo-
% roheptane (16 spins). Different curves correspond the norms of the pro- 
% jection of the density matrix into the subspace of one-, two-, three-,
% etc. spin correlations. The two traces in the lower part of the figure
% correspond to nine- and ten-spin correlations – it is clear that for
% practical simulation purposes, even in the absence of relaxation, only
% correlations of up to eight spins need to be accounted for.
%
% Calculation time: minutes, faster with a GPU.
%
% ilya.kuprov@weizmann.ac.il

function state_spaces_2()
     
% Magnet induction 
sys.magnet=11.7464;

% Isotopes
sys.isotopes={'12C', '12C', '12C', '12C', '12C', '12C', '12C', ...
              '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H', ...
              '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H'};  

% Chemical shifts
inter.zeeman.scalar=cell(1,23);
inter.zeeman.scalar{14}= 1.0092;
inter.zeeman.scalar{15}= 1.0092;
inter.zeeman.scalar{16}= 1.0092;
inter.zeeman.scalar{21}= 1.0092;
inter.zeeman.scalar{22}= 1.0092;
inter.zeeman.scalar{23}= 1.0092;
inter.zeeman.scalar{11}= 4.6834; 
inter.zeeman.scalar{17}= 4.6834;
inter.zeeman.scalar{10}= 0.0000; % actually -184.1865, but does not 
inter.zeeman.scalar{18}= 0.0000; % matter here, and faster when zero
inter.zeeman.scalar{8}=  1.7970;
inter.zeeman.scalar{9}=  1.7970;
inter.zeeman.scalar{13}= 1.6942;
inter.zeeman.scalar{20}= 1.6942;
inter.zeeman.scalar{19}= 1.6370;
inter.zeeman.scalar{12}= 1.6370;

% J-couplings
inter.coupling.scalar=cell(23);  
inter.coupling.scalar{19, 20}= -14.4404; 
inter.coupling.scalar{17, 8}=    2.1802;
inter.coupling.scalar{17, 9}=   10.0564;
inter.coupling.scalar{17, 19}=   4.5308;
inter.coupling.scalar{17, 20}=   7.5739;
inter.coupling.scalar{12, 13}= -14.4404;
inter.coupling.scalar{11, 8}=   10.0564;
inter.coupling.scalar{11, 9}=    2.1802;
inter.coupling.scalar{11, 12}=   4.5308;
inter.coupling.scalar{11, 13}=   7.5739;
inter.coupling.scalar{8,  9}=  -15.1172;
inter.coupling.scalar{10, 8}=   14.0478; 
inter.coupling.scalar{10, 9}=   38.0147; 
inter.coupling.scalar{10, 11}=  49.6307;
inter.coupling.scalar{10, 13}=  18.4317; 
inter.coupling.scalar{10, 12}=  27.7945;
inter.coupling.scalar{18, 9}=   14.0478;
inter.coupling.scalar{18, 8}=   38.0147;
inter.coupling.scalar{18, 17}=  49.6307;
inter.coupling.scalar{18, 20}=  18.4317; 
inter.coupling.scalar{18, 19}=  27.7945; 
inter.coupling.scalar{10, 18}=   1.2295;
inter.coupling.scalar{12,14}=    7.45;
inter.coupling.scalar{12,15}=    7.45;
inter.coupling.scalar{12,16}=    7.45;
inter.coupling.scalar{13,14}=    7.45;
inter.coupling.scalar{13,15}=    7.45;
inter.coupling.scalar{13,16}=    7.45;
inter.coupling.scalar{19,21}=    7.45;
inter.coupling.scalar{19,22}=    7.45;
inter.coupling.scalar{19,23}=    7.45;
inter.coupling.scalar{20,21}=    7.45;
inter.coupling.scalar{20,22}=    7.45;
inter.coupling.scalar{20,23}=    7.45;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=1; bas.manual=false(3,23);
bas.manual(1,[14 15 16 12 13 10 11 8 9])=1;
bas.manual(2,[12 13 10 11 8 9 17 18 19 20])=1;
bas.manual(3,[8 9 17 18 19 20 21 22 23])=1;
bas.sym_group={'S3','S3'};
bas.sym_spins={[14 15 16],[21 22 23]};
bas.longitudinals={'19F'};
bas.projections=1;

% Prevent automatic state dropout
sys.disable={'zte'};

% GPU is useful here
% sys.enable={'gpu'};

% Spinach housekeeping 
spin_system=create(sys,inter); 
spin_system=basis(spin_system,bas);

% Initial condition
rho=state(spin_system,'L+','1H');

% Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Trajectory calculation
traj=evolution(spin_system,H,[],rho,1e-3,1000,'trajectory');

% Trajectory analysis by spin correlation order
kfigure(); trajan(spin_system,traj,'correlation_order'); 
xlim tight; ylim([1e-6 3]); set(gca,'YScale','log');

end

