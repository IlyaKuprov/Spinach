% Commutators of simple operators and superoperators. The test
% calculation is performed three times in the three formalisms
% supported by Spinach.
%
% i.kuprov@soton.ac.uk

function commutation_4()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','235U'};
inter.zeeman.scalar={2.5 1.0};
inter.coupling.scalar{1,2}=10;
inter.coupling.scalar{2,1}=10;

% Preallocate the answer
answer=zeros(3,3,'like',1i);

% Run the tests
formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'};
for n=1:numel(formalisms)
    bas.formalism=formalisms{n};
    bas.approximation='none';
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    CTx=operator(spin_system,'CTx','235U');
    CTy=operator(spin_system,'CTy','235U');
    CTz=operator(spin_system,'CTz','235U');
    CTp=operator(spin_system,'CT+','235U');
    CTm=operator(spin_system,'CT-','235U');
    answer(1,n)=norm(CTz*CTp-CTp*CTz-CTp,'fro');
    answer(2,n)=norm(CTz*CTm-CTm*CTz+CTm,'fro');
    answer(3,n)=norm(CTx*CTy-CTy*CTx-1i*CTz,'fro');
end

% Report the outcome
if norm(answer,'fro')<1e-6
    disp('Cross-formalism commutation test PASSED.');
else
    error('Cross-formalism commutation test FAILED.');
end

end

