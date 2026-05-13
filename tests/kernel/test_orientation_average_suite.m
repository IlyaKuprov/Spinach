% Tests orientation() and average() on small exact cases. Syntax:
%
%                    result=test_orientation_average_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks the zero-Euler-angle orientation contraction against an
% explicit diagonal Wigner sum and verifies that first-order average
% Hamiltonian theory leaves an unmodulated Hamiltonian unchanged.
%
% ilya.kuprov@weizmann.ac.il

function result=test_orientation_average_suite()

% State the rotational-kernel target of the test
result=new_test_result('kernel/orientation_average_suite',...
                       'Orientation and average helpers',...
                       'rotational helper kernels must preserve exact limiting cases.');

% Build a synthetic rank-one rotational basis
Q=cell(1,1);
Q{1}=cell(3,3);
for n=1:3
    for k=1:3
        Q{1}{n,k}=sparse(2,2);
    end
end
Q{1}{1,1}=sparse([1 0;0 0]);
Q{1}{2,2}=sparse([0 1;1 0]);
Q{1}{3,3}=sparse([0 0;0 2]);

% Contract the zero-orientation Hamiltonian
H=orientation(Q,[0 0 0]);
H_ref=Q{1}{1,1}+Q{1}{2,2}+Q{1}{3,3};

% Check the zero-angle Wigner identity path
result=test_close(result,'zero Euler orientation',H,H_ref,1e-14,1e-14,...
                  'wigner(r,0,0,0) is the identity, so only diagonal rotational components remain');

% Build a quiet spin system for average() diagnostics
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Define an unmodulated Hamiltonian decomposition
Hp=sparse(2,2);
H0=sparse([0 1;-1 0]);
Hm=sparse(2,2);
omega=2*pi*1000;

% Run the production average Hamiltonian helper
H_avg=average(spin_system,Hp,H0,Hm,omega,'ah_first_order');

% Check the exact unmodulated limit
result=test_close(result,'unmodulated average',H_avg,H0,1e-14,1e-14,...
                  'with zero positive and negative Fourier components first-order averaging returns H0');

end


