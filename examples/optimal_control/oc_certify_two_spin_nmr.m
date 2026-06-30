% Two-spin magnetic-resonance waveform certification example. A short RF
% waveform is compared with a closed-form reachability bound for a coupled
% Hilbert-space spin pair.
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il

function oc_certify_two_spin_nmr()

% Minimal quiet Spinach object for Hilbert-space propagation
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.sys.disable={};
spin_system.tols.liouv_zero=1e-14;
spin_system.tols.small_matrix=64;
spin_system.tols.dense_matrix=0.5;
spin_system.tols.prop_chop=1e-14;
spin_system.bas.formalism='zeeman-hilb';

% Spin-half operators
S=pauli(2);
E=speye(2);

% Build two-spin product operators
Lx=kron(S.x,E)+kron(E,S.x);
Ly=kron(S.y,E)+kron(E,S.y);
Lz=kron(S.z,E)+kron(E,S.z);
Lzz=kron(S.z,S.z);

% Optimal-control metadata
spin_system.control.operators={Lx,Ly};
spin_system.control.rho_init={Lz};
spin_system.control.rho_targ={Lx};
spin_system.control.pwr_levels=2*pi*1000;
spin_system.control.pulse_dt=20e-6*ones(1,4);
spin_system.control.drifts={{2*pi*12*Lzz}};
spin_system.control.integrator='rectangle';
spin_system.control.dead_time=0;
spin_system.control.prefix=[];
spin_system.control.suffix=[];
spin_system.control.keyholes=cell(1,numel(spin_system.control.pulse_dt));

% Candidate normalised RF waveform and certificate
problem.waveform=[0.35 0.65 0.65 0.35; 0.00 0.20 -0.20 0.00];
certificate=oc_certify(spin_system,problem);

% Report the result
disp(certificate);

end
