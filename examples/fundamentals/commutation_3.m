% Commutators of simple operators and superoperators. All output 
% should be close to zero. The same calculation is performed four 
% times in the four formalisms supported by Spinach.
%
% The calculation should return an 8x3 zero matrix.
%
% i.kuprov@soton.ac.uk

function commutation_3()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','235U'};
inter.zeeman.scalar={2.5 1.0};
inter.coupling.scalar{1,2}=10;
inter.coupling.scalar{2,1}=10;

% Preallocate the answer
answer=zeros(7,3,'like',1i);

% Run the tests
formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'};
for n=1:numel(formalisms)
    bas.formalism=formalisms{n};
    bas.approximation='none';
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    Up=operator(spin_system,'L+','235U');
    Um=operator(spin_system,'L-','235U');
    Uz=operator(spin_system,'Lz','235U');
    HpUp=operator(spin_system,{'L+','L+'},{1,2});
    HmUm=operator(spin_system,{'L-','L-'},{1,2});
    Hz=operator(spin_system,'Lz','1H');
    Ux=(Up+Um)/2; Uy=(Up-Um)/2i;
    answer(1,n)=norm(Uz*Up-Up*Uz-Up,'fro');
    answer(2,n)=norm(Uz*Um-Um*Uz+Um,'fro');
    answer(3,n)=norm(Ux*Uy-Uy*Ux-1i*Uz,'fro');
    answer(4,n)=norm(Uz*HpUp-HpUp*Uz-HpUp,'fro');
    answer(5,n)=norm(Hz*HpUp-HpUp*Hz-HpUp,'fro');
    answer(6,n)=norm(Uz*HmUm-HmUm*Uz+HmUm,'fro');
    answer(7,n)=norm(Hz*HmUm-HmUm*Hz+HmUm,'fro');
end

% Display the answers
disp(abs(answer)>1e-6);

end

