% Relaxation superoperators extracted from stochastic Liouville equation
% generators by Markovian back-projection onto the spin Liouville space
% (sle2rlx.m), compared with Redfield theory (relaxation.m) for the same
% spin systems and the same rotational correlation times: the nitroxide
% radical of sle_esr_nitroxide_1.m and the 15N-1H pair with dipolar and
% CSA interactions of sle_nmr_dd_csa.m. Both calculations are done in the
% laboratory frame with the dynamic frequency shifts retained. Pertinent
% matrix elements of the two superoperators are printed to the console;
% the residual differences are the higher-order corrections to Redfield
% theory that the full SLE memory kernel contains.
%
% Calculation time: about a minute
%
% talos@spindynamics.org

function sle2rlx_vs_redfield()

% Nitroxide radical
options.no_xyz=1;
[sys,inter]=g2spinach(gparse('../standard_systems/nitroxide.log'),...
                           {{'E','E'},{'N','14N'}},[0 0],options);
sys.magnet=3.5;
sys.tols.prox_cutoff=4.0;

% Redfield theory, lab frame, dynamic frequency shifts kept
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.rlx_dfs='keep';
inter.tau_c={5e-11};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation={'none'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Redfield superoperator
R_rf=relaxation(spin_system);

% SLE generator in the lab frame
parameters.spins={'E'};
parameters.max_rank=4;
parameters.tau_c=5e-11;
capture=@(spin_system,parameters,H,R,K)H+1i*R+1i*K;
L=gridfree(spin_system,capture,parameters,'labframe');

% Lattice equilibrium is the rank zero Wigner function
P=sparse(1,1,1,size(L,1)/size(R_rf,1),1);

% Back-projection onto the spin Liouville space
R_sle=sle2rlx(spin_system,L,P,P);

% Pertinent matrix elements
disp('Nitroxide radical, rad/s:');
Ez=state(spin_system,'Lz','E'); Ep=state(spin_system,'L+','E');
EpNz=state(spin_system,{'L+','Lz'},{1,2});
disp_element('<Ez|R|Ez>/<Ez|Ez>',Ez,Ez,R_rf,R_sle);
disp_element('<E+|R|E+>/<E+|E+>',Ep,Ep,R_rf,R_sle);
disp_element('<E+Nz|R|E+>/|E+Nz||E+|',EpNz,Ep,R_rf,R_sle);
disp(['||R_sle-R_rf||/||R_rf|| = ' num2str(norm(R_sle-R_rf,'fro')/norm(R_rf,'fro'))]);

% 15N-1H spin pair
clear('sys','inter','parameters');
sys.magnet=14.1;
sys.isotopes={'15N','1H'};
inter.zeeman.matrix{1}=[14.89    34.35     0.00
                        34.35   145.74     0.00
                         0.00     0.00   150.51];
inter.zeeman.matrix{2}=[30.75    1.31    0.00
                         1.31   22.65    0.00
                         0.00    0.00   11.80];
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=10;
inter.coordinates={[-0.451455   -0.678015    0.000000]
                   [-1.475290   -0.641823    0.000000]};
sys.tols.prox_cutoff=4.0;

% Redfield theory, lab frame, dynamic frequency shifts kept
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.rlx_dfs='keep';
inter.tau_c={5e-9};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Redfield superoperator
R_rf=relaxation(spin_system);

% SLE generator in the lab frame
parameters.spins={'15N'};
parameters.max_rank=4;
parameters.tau_c=5e-9;
L=gridfree(spin_system,capture,parameters,'labframe');

% Lattice equilibrium is the rank zero Wigner function
P=sparse(1,1,1,size(L,1)/size(R_rf,1),1);

% Back-projection onto the spin Liouville space
R_sle=sle2rlx(spin_system,L,P,P);

% Pertinent matrix elements
disp('15N-1H spin pair, rad/s:');
Nz=state(spin_system,'Lz','15N'); Np=state(spin_system,'L+','15N');
Hz=state(spin_system,'Lz','1H');  NpHz=state(spin_system,{'L+','Lz'},{1,2});
disp_element('<Nz|R|Nz>/<Nz|Nz>',Nz,Nz,R_rf,R_sle);
disp_element('<N+|R|N+>/<N+|N+>',Np,Np,R_rf,R_sle);
disp_element('<Hz|R|Nz>/|Hz||Nz|',Hz,Nz,R_rf,R_sle);
disp_element('<N+Hz|R|N+>/|N+Hz||N+|',NpHz,Np,R_rf,R_sle);
disp(['||R_sle-R_rf||/||R_rf|| = ' num2str(norm(R_sle-R_rf,'fro')/norm(R_rf,'fro'))]);

end

% Prints a normalised matrix element of both superoperators
function disp_element(label,bra,ket,R_rf,R_sle)
el_rf=(bra'*R_rf*ket)/(norm(bra)*norm(ket));
el_sle=(bra'*R_sle*ket)/(norm(bra)*norm(ket));
fprintf('%-26s Redfield %+.6e%+.6ei  SLE %+.6e%+.6ei\n',...
        label,real(el_rf),imag(el_rf),real(el_sle),imag(el_sle));
end

