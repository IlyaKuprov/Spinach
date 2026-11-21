% Basal kh-divacancy free induction at the zero-field clock transition.
% Reconstructs the seed-1 nuclear configuration used for
% Fig. 4(c) of Onizhuk and Galli, doi:10.1002/adts.202100254, using
% mixed nuclear states and 40 Angstrom bath and 8 Angstrom pair cutoffs.
%
% Returns the normalised complex gCCE2 coherence and diagnostics,
% including the gCCE1 coherence and numerical error bounds.
%
% Calculation time: about 3 minutes on one Xeon Platinum 8270 worker.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=gcce_sic_clock()

% Read the fixed atomistic configuration and tensors in Angstrom and Hz
bath=load(fullfile(fileparts(mfilename('fullpath')),'sic_clock_kh_bath.mat'));
nsites=numel(bath.isotopes)+1;
sys.magnet=0;
sys.isotopes=[{'E3'} bath.isotopes];
sys.output='hush';
sys.parallel={'processes',1};
inter.coupling.matrix=cell(nsites);
inter.zeeman.matrix=cell(1,nsites);

% Spin-one axial and transverse zero-field splitting in Hz
axial=1.334e9;
transverse=18.4e6;
inter.coupling.matrix{1,1}=diag([-axial/3+transverse -axial/3-transverse 2*axial/3]);

% Full hyperfine interactions
for n=2:nsites

    % PyCCE stores nuclear index before electronic index
    inter.coupling.matrix{1,n}=bath.hyperfine(:,:,n-1).';

end

% Nuclear dipolar interactions
for n=1:size(bath.edges,1)
    inter.coupling.matrix{bath.edges(n,1),bath.edges(n,2)}=bath.dipolar(:,:,n);
end

% Prepare and detect the coherence between the mS=0 and even defect eigenstates
ground=[0;1;0];
excited=[1;0;1]/sqrt(2);
prepared=(ground+excited)/sqrt(2);
parameters.observed=1;
parameters.rho0=prepared*prepared';
parameters.coil=ground*excited';
parameters.normalized=true;
parameters.times=linspace(0,1e-3,501);

% Mixed-state generalised expansion with the 8 Angstrom pair graph
adjacency=sparse(bath.edges(:,1)-1,bath.edges(:,2)-1,true,nsites-1,nsites-1);
bas.formalism='zeeman-hilb';
parameters.cce.approximation='generalized';
parameters.cce.basis=bas;
parameters.cce.order=2;
parameters.cce.graph=struct('mode','explicit','adjacency',adjacency|adjacency.',...
                            'approximate',true);
parameters.cce.bath_state=struct('mode','mixed_product');
parameters.cce.meanfield='none';

% Absolute budget covers the measured worst-case order-two roundoff bound
parameters.cce.numerics=struct('atol',0.03,'rtol',1e-8,...
                               'max_memory_bytes',2^31);

% Calculate the signal without building the full bath basis
spin_system=create(sys,inter);
[answer,diagnostics]=clexpn(spin_system,@free_evolution,parameters,'labframe');

% Compare the first two cluster orders
kfigure();
plot(1e3*parameters.times,abs(squeeze(diagnostics.reduction.by_order(:,1,2:3))));
kxlabel('time / ms');
kylabel('normalised coherence magnitude');
klegend({'gCCE1','gCCE2'});
kgrid();

end


