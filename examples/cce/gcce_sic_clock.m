% Basal kh-divacancy free induction at the zero-field clock transition.
% Reconstructs the seed-1 nuclear configuration used for
% Fig. 4(c) of Onizhuk and Galli, doi:10.1002/adts.202100254, using
% mixed nuclear states and 40 Angstrom bath and 8 Angstrom pair cutoffs.
%
% Calculation time: about 3 minutes on one Xeon Platinum 8270 worker.
%
% ilya.kuprov@weizmann.ac.il

function gcce_sic_clock()

% Read the fixed atomistic configuration and tensors in Angstrom and Hz
folder=fileparts(mfilename('fullpath'));
sites=readtable(fullfile(folder,'sic_clock_kh_sites.csv'),'TextType','string');
pairs=readtable(fullfile(folder,'sic_clock_kh_pairs.csv'));
hyperfine=reshape(sites{:,6:end}.',3,3,[]);
dipolar=reshape(pairs{:,3:end}.',3,3,[]);
nsites=height(sites)+1;
sys.magnet=0;
sys.isotopes={'E3'};
sys.isotopes(sites.spin)=cellstr(sites.isotope);
sys.output='hush';
inter.coupling.matrix=cell(nsites);
inter.zeeman.matrix=cell(1,nsites);

% Spin-one axial and transverse zero-field splitting in Hz
axial=1.334e9;
transverse=18.4e6;
inter.coupling.matrix{1,1}=diag([-axial/3+transverse -axial/3-transverse 2*axial/3]);

% Full hyperfine interactions
for n=1:height(sites)

    % PyCCE stores nuclear index before electronic index
    inter.coupling.matrix{1,sites.spin(n)}=hyperfine(:,:,n).';

end

% Nuclear dipolar interactions
for n=1:height(pairs)
    inter.coupling.matrix{pairs.spin_1(n),pairs.spin_2(n)}=dipolar(:,:,n);
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
adjacency=sparse(pairs.spin_1-1,pairs.spin_2-1,true,nsites-1,nsites-1);
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
[~,details]=clexpn(spin_system,@free_evolution,parameters,'labframe');

% Compare the first two cluster orders
kfigure();
plot(1e3*parameters.times,abs(squeeze(details.reduction.by_order(:,1,2:3))));
kxlabel('time / ms');
kylabel('normalised coherence magnitude');
klegend({'gCCE1','gCCE2'});
kgrid();

end


