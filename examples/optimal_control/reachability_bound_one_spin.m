% One-spin reachability certificate example. A bounded y-pulse attempts
% to rotate Lx into Lz; reachable_bound() gives a rigorous upper bound
% on the target fidelity at the requested duration.
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il

function reachability_bound_one_spin()

% Spin-half operators
S=pauli(2);

% Define a two-level control problem
drift=sparse(2,2);
controls={S.y};
rho_init=S.x;
rho_targ=S.z;
amplitude_bound=100;
duration=0.01;

% Certify the requested duration
certificate=reachable_bound(drift,controls,rho_init,rho_targ,...
                            amplitude_bound,duration);

% Report the result
disp(certificate);

end

