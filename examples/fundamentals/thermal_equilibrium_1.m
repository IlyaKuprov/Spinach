% Observables at thermal equilibrium using the three formalisms 
% supported by Spinach kernel, tested against known answers.
%
% i.kuprov@soton.ac.uk

function thermal_equilibrium_1()

% Spin system parameters
sys.magnet=14.1;
sys.isotopes={'E','1H','14N','15N'};
inter.zeeman.scalar={2.0023 1.0 2.0 3.0};
inter.coupling.scalar=cell(length(sys.isotopes));
inter.coupling.scalar{1,2}=1e6;
inter.coupling.scalar{2,3}=1e6;
inter.coupling.scalar{1,3}=1e3;
inter.coupling.scalar{3,4}=1e3;
inter.coupling.scalar{1,4}=1e2;
inter.temperature=4.2;

% Preallocate the answers
eq_mags_spinach=zeros(4,3,'like',1i);
eq_mags_textbook=zeros(4,1,'like',1i);

% Get numerical equilibrium magnetisation
formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'};
for n=1:numel(formalisms)
    bas.formalism=formalisms{n};
    bas.approximation='none';
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    rho=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'));
    eq_mags_spinach(1,n)=trace(state(spin_system,'Lz','E')'*rho);
    eq_mags_spinach(2,n)=trace(state(spin_system,'Lz','1H')'*rho);
    eq_mags_spinach(3,n)=trace(state(spin_system,'Lz','14N')'*rho);
    eq_mags_spinach(4,n)=trace(state(spin_system,'Lz','15N')'*rho);
end

% Get analytical equilibrium magnetisation
[~,P]=levelpop('E',sys.magnet,inter.temperature);
eq_mags_textbook(1)=0.5*P(1)-0.5*P(2);
[~,P]=levelpop('1H',sys.magnet,inter.temperature);
eq_mags_textbook(2)=0.5*P(1)-0.5*P(2);
[~,P]=levelpop('14N',sys.magnet,inter.temperature);
eq_mags_textbook(3)=P(1)-P(3);
[~,P]=levelpop('15N',sys.magnet,inter.temperature);
eq_mags_textbook(4)=0.5*P(1)-0.5*P(2);

% Display the answers
if any(abs(eq_mags_spinach-eq_mags_textbook)>1e-5,'all')
    error('Cross-formalism thermal equilibrium test FAILED.');
else
    disp('Cross-formalism thermal equilibrium test PASSED.');
end

end

