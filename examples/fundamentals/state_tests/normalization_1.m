% Internal consistency test for the state vectors and matrices
% in each of the three formalisms supported by Spinach.
%
% i.kuprov@soton.ac.uk

function normalization_1()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','235U'};
inter.zeeman.scalar={2.5 1.0};
inter.coupling.scalar{1,2}=10;
inter.coupling.scalar{2,1}=10;
spin_system=create(sys,inter);

% Preallocate the answer
norm_diffs=zeros(6,3,'like',1i);

% Compute norm differences
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
    norm_diffs(1,n)=norm(full(Ux))-norm(full(Uy));
    norm_diffs(2,n)=norm(full(Uy))-norm(full(Uz));
    norm_diffs(3,n)=norm(full(Uz))-norm(full(Ux));
    norm_diffs(4,n)=norm(full(Hx))-norm(full(Hy));
    norm_diffs(5,n)=norm(full(Hy))-norm(full(Hz));
    norm_diffs(6,n)=norm(full(Hz))-norm(full(Hx));
end

% Display the answers
if any(abs(norm_diffs)>1e-6,'all')
    error('Internal norm consistency test FAILED.');
else
    disp('Internal norm consistency test PASSED.');
end

end

