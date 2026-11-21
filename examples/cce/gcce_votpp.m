% Field-dependent energy levels of the electron–vanadium-51 register
% in VO(TPP), using Eq. (10) and Table I to reconstruct Figure 2 of
% doi:10.1063/5.0278659. For Figure 2 the field is transverse to the unique axis
% of the g, hyperfine, and quadrupolar tensors. The proton Hahn echo
% (middle of Figure 3) uses projected CCE-2 with the source crystal
% coordinates and tensors, retaining the entire proton bath inside 20 A.
% For Figure 3 the companion input swaps x and z in g and A only;
% the scalar vanadium quadrupolar splitting remains axial along z.
% Syntax:
%
%                    [energies,fields,echo,times]=gcce_votpp()
%
% Rows of energies contain the 16 ascending energy levels in GHz
% at the corresponding entries of the field row vector fields in tesla.
% echo and times are rows of complex normalised coherence and total
% Hahn echo times (seconds) at 0.33 T.
%
% ilya.kuprov@weizmann.ac.il

function [energies,fields,echo,times]=gcce_votpp()

% The unique axis is x in the crystal frame; the applied field is along z
sys.isotopes={'E','51V'};
sys.output='hush';
sys.parallel={'processes',1};
inter.zeeman.matrix={diag([1.968 1.984 1.984]),zeros(3)};
inter.coupling.matrix=cell(2);
inter.coupling.matrix{1,2}=diag([-473 -166 -166])*1e6;
inter.coupling.matrix{2,2}=-0.35e6*diag([2 -1 -1])/3;

% Compile the observed-only gCCE Hamiltonian at the field endpoints
field_ends=[0 0.15];
parameters.observed=[1 2];
parameters.times=0;
parameters.cce.approximation='generalized';
parameters.cce.order=0;
parameters.cce.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^20);
H=cell(1,2);
for k=1:2

    % Convert physical tensors and select the empty bath cluster
    sys.magnet=field_ends(k);
    spin_system=create(sys,inter);
    model=cce_compile(spin_system,parameters,'labframe');
    H{k}=cce_hamiltonian(model,[],struct());

end

% Interpolate the Hamiltonian, not its eigenvalues, because Zeeman terms are linear in field
zeeman=(H{2}-H{1})/field_ends(2);
fields=linspace(field_ends(1),field_ends(2),100);
energies=zeros(numel(fields),size(H{1},1));
for k=1:numel(fields)

    % Sort the energy levels in spectroscopic units
    energies(k,:)=sort(eig(full(H{1}+fields(k)*zeeman))).'/(2*pi*1e9);

end

% Plot the field dependence
kfigure();
plot(fields,energies);
kxlabel('magnetic field / T');
kylabel('energy / GHz');
kgrid();

% Calculate and display the proton-limited molecular echo
[echo,times]=votpp_proton_echo();
kfigure();
plot(1e6*times,real(echo));
kxlabel('total echo time / microseconds');
kylabel('normalised coherence');
kgrid();

end


