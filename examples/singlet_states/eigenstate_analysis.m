% Stationary state analysis for the spin system of allyl pyruvate,
% finding out which component of the singlet state commutes with
% the drift Hamiltonian.
%
% i.kuprov@soton.ac.uk

function eigenstate_analysis()

% Get the spin system from Anu's fits
[sys,inter]=allyl_pyruvate({'1H','13C'});

% Set the magnet
sys.magnet=14.1;

% Spinach housekeeping
spin_system=create(sys,inter);

% Pick out the required 13C isotopomer
subsystems=dilute(spin_system,'13C',1);
spin_system=subsystems{4};

% Generate the basis
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=basis(spin_system,bas);

% Get isotropic Hamiltonian
spin_system=assume(spin_system,'nmr');
H=hamiltonian(spin_system); H=full(H);

% Tidy up rounding errors
H=(H+ctranspose(H))/2;

% Get the singlet state
rho=singlet(spin_system,3,4);

% Report the norm
report(spin_system,['Singlet state norm: ' ...
                    num2str(norm(full(rho),2))]);

% Remove the part that does not commute with H
[EvH,~]=eig(H); rho_inv=remncomm(rho,EvH);

% Remove the unit part
rho_inv=remtrace(rho_inv);

% Report the norm
report(spin_system,['...of which commutes with H: ' ...
                     num2str(norm(full(rho_inv),2))]);

% Project out the ZZ component
zz_state=state(spin_system,{'Lz','Lz'},{3 4});
zz_state=zz_state/norm(zz_state,'fro');
rho_inv=rho_inv-zz_state*trace(zz_state'*rho_inv);

% Report the norm
report(spin_system,['...of which is not ZZ: ' ...
                     num2str(norm(full(rho_inv),2))]);

%% Apply offset and spin-lock, repeat process 

% Hamiltonian offset
parameters.offset=2850;
parameters.spins={'1H'};
H=frqoffset(spin_system,H,parameters);

% Spin-lock term at 1 kHz
H=H+2*pi*1000*operator(spin_system,'Lx','1H');

% Get the singlet state
rho=singlet(spin_system,3,4);

% Report the norm
report(spin_system,['Singlet state norm: ' ...
                    num2str(norm(full(rho),2))]);

% Tidy up rounding errors
H=(H+ctranspose(H))/2;

% Remove the part that does not commute with H
[EvH,~]=eig(H); rho_inv=remncomm(rho,EvH);

% Remove the unit part
rho_inv=remtrace(rho_inv);

% Report the norm
report(spin_system,['...of which commutes with H: ' ...
                     num2str(norm(full(rho_inv),2))]);

% Project out the ZZ component
zz_state=state(spin_system,{'Lz','Lz'},{3 4});
zz_state=zz_state/norm(zz_state,'fro');
rho_inv=rho_inv-zz_state*trace(zz_state'*rho_inv);

% Report the norm
report(spin_system,['...of which is not ZZ: ' ...
                     num2str(norm(full(rho_inv),2))]);

end

