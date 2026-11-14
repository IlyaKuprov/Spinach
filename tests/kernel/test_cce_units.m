% Tests CCE canonical units against independent spin-half formulas. Syntax:
%
%                         result=test_cce_units()
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_units()

% State the independent numerical targets
result=new_test_result('kernel/cce_units','CCE canonical units',...
                       'Check signed Zeeman, scalar, dipolar and ordered-cell units.');

% Define explicit model choices without invoking an interaction conversion twice
parameters.central_spins=1; parameters.bath_spins=2;
parameters.assumptions='labframe'; parameters.orientation=[0 0 0];
parameters.interaction_cutoff=Inf;
parameters.numerics=struct('atol',1e-7,'rtol',1e-12,'max_memory_bytes',1e8);
sys.isotopes={'1H','13C'}; sys.magnet=0; sys.output='hush';
inter.zeeman.scalar={0,0}; inter.coupling.scalar={0,137;0,0};

% Compare scalar Hz conversion with the explicit singlet/triplet matrix
model=cce_prepare(create(sys,inter),parameters);
[H,meta]=cce_hamiltonian(model,2,struct());
reference=2*pi*137*[1/4 0 0 0;0 -1/4 1/2 0;0 1/2 -1/4 0;0 0 0 1/4];
result=test_close(result,'scalar Hz to rad/s',...
                  H,reference,meta.roundoff_bound,0,'Independent spin-half scalar matrix.');

% Compare coordinate dipoles with the SI point-dipole expression
inter=struct('zeeman',struct('scalar',{{0,0}}));
inter.coordinates={[0 0 0],[0 0 2]};
model=cce_prepare(create(sys,inter),parameters);
[H,meta]=cce_hamiltonian(model,2,struct());
gamma_carbon=0.70236944*5.0507837393e-27/(0.5*1.054571817e-34);
dipole=1e-7*1.054571817e-34*267522187.08*gamma_carbon/(2e-10)^3;
reference=dipole*[-1/2 0 0 0;0 1/2 1/2 0;0 1/2 1/2 0;0 0 0 -1/2];
result=test_close(result,'dipole SI units',...
                  H,reference,0,2e-7,'Independent gamma/hbar formula and Angstrom conversion.');

% Check one full ordered cell versus two transposed half cells
system=model.spin_system; tensor=system.inter.coupling.matrix{1,2}+...
                                  system.inter.coupling.matrix{2,1}.';
system.inter.coupling.matrix{1,2}=tensor;
system.inter.coupling.matrix{2,1}=zeros(3);
one=cce_hamiltonian(cce_prepare(system,parameters),2,struct());
result=test_close(result,'ordered half cells',...
                  H,one,meta.roundoff_bound,0,'Coordinate tensors must not double-count.');

% Check the nuclear Zeeman sign and a nonzero electron field independently
for isotope={'1H','E'}
    sys.isotopes=isotope; sys.magnet=0.07;
    parameters.bath_spins=[];
    model=cce_prepare(create(sys,struct()),parameters);
    [H,meta]=cce_hamiltonian(model,[],struct());
    if strcmp(isotope{1},'1H'), gamma=267.52218744e6; else, gamma=-1.76085963023e11; end
    reference=-gamma*sys.magnet*diag([1/2 -1/2]);
    result=test_close(result,[isotope{1} ' signed Zeeman'],...
                  H,reference,0,2e-7,'Independent signed gyromagnetic ratio.');
    result=test_true(result,'sign mutation rejected',...
                  norm(H+reference,'fro')>100*meta.roundoff_bound,'The opposite Zeeman sign is distinguishable.');
end

% Reject unknown physical choices rather than silently accepting a typo
parameters.unrecognised=1; rejected=false;
try
    cce_prepare(model.spin_system,parameters);
catch
    rejected=true;
end
result=test_true(result,'unknown field rejected',...
                  rejected,'Model choices must be explicit.');

end


