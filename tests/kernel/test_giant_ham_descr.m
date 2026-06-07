% Tests the giant spin Hamiltonian descriptor route. Syntax:
%
%                    result=test_giant_ham_descr()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that high-rank giant spin Hamiltonian terms assembled
% through the descriptor route match direct spherical-tensor assembly.
%
% ilya.kuprov@weizmann.ac.il

function result=test_giant_ham_descr()

% Announce the test target
fprintf('TESTING: Giant spin Hamiltonian descriptor\n');

% State the Hamiltonian target of the test
result=new_test_result('kernel/giant_ham_descr',...
                       'Giant spin Hamiltonian descriptor',...
                       'high-rank giant spin terms must match direct spherical-tensor assembly.');

% Check complete giant spin terms
result=local_case(result,'labframe','strong',[0.41 0.29 0.13]);

% Check secular giant spin terms
result=local_case(result,'deer-zz','secular',[0.17 0.39 0.51]);

end


% Checks one giant spin Hamiltonian assumption
function result=local_case(result,assumption,strength,euler_angles)

% Build a compact high-rank giant spin system
sys.magnet=0;
sys.isotopes={'E8'};
inter.zeeman.scalar={0};
inter.giant.coeff={{[0.07e6 -0.01e6 0.03e6],...
                    [0.02e6 0.01e6 -0.03e6 0.04e6 -0.02e6],...
                    [0.05e6 -0.02e6 0.03e6 0.04e6 -0.01e6 0.02e6 -0.03e6]}};
inter.giant.euler={{[0.11 0.13 0.17],...
                    [0.19 0.23 0.29],...
                    [0.31 0.37 0.41]}};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Apply the requested Hamiltonian assumption
spin_system=assume(spin_system,assumption);

% Build the production Hamiltonian
[I,Q]=hamiltonian(spin_system);
H_obs=I+orientation(Q,euler_angles);

% Build the direct reference Hamiltonian
H_ref=local_ref(spin_system,euler_angles);

% Check that high-rank components are present
result=test_close(result,['giant spin rank depth ' strength],...
                  double(numel(Q)>=3),1,0,0,...
                  'rank-three giant spin terms must extend the rotational basis');

% Check the oriented Hamiltonian
result=test_close(result,['giant spin Hamiltonian ' strength],H_obs,H_ref,...
                  1e-7,1e-12,...
                  'descriptor assembly must match direct spherical-tensor assembly');

end


% Builds a direct spherical-tensor giant spin Hamiltonian reference
function H_ref=local_ref(spin_system,euler_angles)

% Start from a zero Hamiltonian
H_ref=mprealloc(spin_system,0);

% Loop over spins
for n=1:spin_system.comp.nspins

    % Loop over spherical ranks
    for r=1:numel(spin_system.inter.giant.coeff{n})

        % Compute the Wigner matrix
        W=wigner(r,euler_angles(1),euler_angles(2),euler_angles(3));

        % Process the giant spin assumption
        switch spin_system.inter.giant.strength{n}

            case 'strong'

                % Loop over spherical tensor projections
                for k=1:(2*r+1)

                    % Contract coefficient components
                    coeff=W(k,:)*spin_system.inter.giant.coeff{n}{r}(:);

                    % Add the operator contribution
                    if abs(coeff)>spin_system.tols.liouv_zero
                        ist_spec=['T' num2str(r) ',' num2str(r-k+1)];
                        H_ref=H_ref+coeff*operator(spin_system,{ist_spec},{n});
                    end

                end

            case 'secular'

                % Contract coefficient components
                coeff=W(r+1,:)*spin_system.inter.giant.coeff{n}{r}(:);

                % Add the operator contribution
                if abs(coeff)>spin_system.tols.liouv_zero
                    ist_spec=['T' num2str(r) ',0'];
                    H_ref=H_ref+coeff*operator(spin_system,{ist_spec},{n});
                end

            case 'ignore'

                % Move to the next term
                continue

        end

    end

end

% Match orientation.m Hermitisation
H_ref=(H_ref+H_ref')/2;

end


