% Two-spin magnetic-resonance reachability certificate. A bounded RF field
% attempts to rotate longitudinal magnetisation into transverse
% magnetisation; reachable_bound() gives a rigorous conservative upper
% bound on the target overlap at the requested duration.
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il

function reachability_bound_two_spin_nmr()

% Spin-half operators
S=pauli(2);
E=speye(2);

% Build two-spin product operators
Lx=kron(S.x,E)+kron(E,S.x);
Ly=kron(S.y,E)+kron(E,S.y);
Lz=kron(S.z,E)+kron(E,S.z);
Lzz=kron(S.z,S.z);

% Define a simple coupled-spin Hamiltonian
drift=2*pi*12*Lzz;
controls={Lx,Ly};
rho_init=Lz;
rho_targ=Lx;
amplitude_bound=2*pi*1000;
duration=80e-6;

% Certify the requested duration
certificate=reachable_bound(drift,controls,rho_init,rho_targ,...
                            amplitude_bound,duration);

% Report the result
disp(certificate);

end
