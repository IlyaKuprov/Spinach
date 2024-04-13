% Thermal equilibrium states, using all the different formalisms 
% supported by Spinach kernel.
%
% i.kuprov@soton.ac.uk

function thermal_equilibrium_2()

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

% Get thermal equilibrium states at finite temperature
formalisms={'zeeman-hilb','zeeman-liouv','sphten-liouv'};
rho=cell(1,numel(formalisms));
for n=1:numel(formalisms)
    bas.formalism=formalisms{n}; bas.approximation='none';
    spin_system=create(sys,inter); spin_system=basis(spin_system,bas);
    rho{n}=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'));
end

% Zeeman-liouv to zeeman-hilb
rho{2}=reshape(rho{2},[24 24]);

% Sphten-liouv to zeeman-hilb
mult_factor=sqrt(prod(spin_system.comp.mults));
rho{3}=sphten2zeeman(spin_system)*rho{3};
rho{3}=reshape(rho{3},[24 24])/mult_factor;

% Check the differences
if (norm(rho{1}-rho{2},2)>1e-8)||(norm(rho{2}-rho{3},2)>1e-8)
    error('Cross-formalism thermal equilibrium test FAILED.');
else
    disp('Cross-formalism thermal equilibrium test PASSED.');
end

end

