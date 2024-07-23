% Commutators of simple operators and superoperators. The test
% calculation is performed three times in the three formalisms
% supported by Spinach.
%
% i.kuprov@soton.ac.uk

function commutation_2()

% System specification
sys.magnet=0;
sys.isotopes={'235U'};
inter.zeeman.scalar={0};

% Preallocate the answer
answer=zeros(3,3,'like',1i);

% Run the tests
formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'};
for n=1:numel(formalisms)
    bas.formalism=formalisms{n};
    bas.approximation='none';
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    Lp=operator(spin_system,'L+','235U');
    Lm=operator(spin_system,'L-','235U');
    Lz=operator(spin_system,'Lz','235U');
    Lx=(Lp+Lm)/2; Ly=(Lp-Lm)/2i;
    answer(1,n)=norm(Lz*Lp-Lp*Lz-Lp,'fro');
    answer(2,n)=norm(Lz*Lm-Lm*Lz+Lm,'fro');
    answer(3,n)=norm(Lx*Ly-Ly*Lx-1i*Lz,'fro');
end

% Report the outcome
if norm(answer,'fro')<1e-6
    disp('Cross-formalism commutation test PASSED.');
else
    error('Cross-formalism commutation test FAILED.');
end

end

