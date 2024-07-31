% Internal consistency test for the state vectors and matrices.
% Checks that the inner products are consistent between the for-
% malisms supported by Spinach. Output should be a 6x3 matrix
% with identical columns.
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de

function normalization_2()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','235U'};
inter.zeeman.scalar={2.5 1.0};
inter.coupling.scalar{1,2}=10;
inter.coupling.scalar{2,1}=10;
spin_system=create(sys,inter);

% Preallocate the answer
norms=zeros(6,3,'like',1i);

% Get the norms
formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'};
for n=1:numel(formalisms)
    bas.formalism=formalisms{n};
    bas.approximation='none';
    spin_system=basis(spin_system,bas);
    Ux=state(spin_system,'Lx','235U');
    Uy=state(spin_system,'Ly','235U');
    Uz=state(spin_system,'Lz','235U');
    Hx=state(spin_system,'Lx','1H');
    Hy=state(spin_system,'Ly','1H');
    Hz=state(spin_system,'Lz','1H');
    norms(1,n)=trace(full(Ux)'*full(Ux));
    norms(2,n)=trace(full(Uy)'*full(Uy));
    norms(3,n)=trace(full(Uz)'*full(Uz));
    norms(4,n)=trace(full(Hx)'*full(Hx));
    norms(5,n)=trace(full(Hy)'*full(Hy));
    norms(6,n)=trace(full(Hz)'*full(Hz));
end

% Run the tests
mult_factor=prod(spin_system.comp.mults);
if norm(norms(:,1)-norms(:,2),1)>1e-6||...
   norm(norms(:,1)-mult_factor*norms(:,3),1)>1e-6
    error('Cross-formalism state norm test FAILED.');
else
    disp('Cross-formalism state norm test PASSED.');
end

end

