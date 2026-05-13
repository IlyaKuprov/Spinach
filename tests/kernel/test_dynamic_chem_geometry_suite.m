% Tests deterministic chemistry and geometry utility helpers. Syntax:
%
%                    result=test_dynamic_chem_geometry_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks lattice construction, geometric measurements, coupling
% extraction, chemical shifts, nearest-neighbour lookup, and tensor helpers.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_chem_geometry_suite()

% State the chemistry and geometry target of the test
result=new_test_result('kernel/dynamic_chem_geometry_suite',...
                       'Chemistry and geometry utilities',...
                       'small chemistry and geometry helpers must preserve exact coordinate, tensor, and metadata formulae.');

% Build a small spin-system descriptor for metadata helpers
spin_system=local_geometry_spin_system();

% Check simple cubic lattice construction and periodic vectors
[sys,inter]=cubic_lattice('13C',2,2);
result=test_true(result,'cubic_lattice isotope count',numel(sys.isotopes)==8&&all(strcmp(sys.isotopes,'13C')),...
                 'a 2x2x2 cubic lattice must produce eight identical isotope entries');
result=test_close(result,'cubic_lattice coordinate order',inter.coordinates{2},[0 0 2],1e-15,1e-15,...
                  'sub2ind ordering places the second atom one spacing along the z coordinate');
result=test_true(result,'cubic_lattice pbc vectors',isequal(inter.pbc,{[4 0 0],[0 4 0],[0 0 4]}),...
                 'periodic boundary vectors must be spacing times the number of periods');

% Check a signed right-angle dihedral from four Cartesian points
phi=dihedral([1 0 0],[0 0 0],[0 1 0],[0 1 1]);
result=test_close(result,'dihedral right angle',phi,-90,1e-12,1e-12,...
                  'the selected A-B-C-D geometry has a minus ninety-degree torsion in Spinach convention');

% Check Cartesian point-cloud binning on a regular grid
coords=[0.25 0.25 0.25;0.75 0.25 0.25;1.5 0.5 0.5];
density=xyz2pd(coords,[0 1],[0 1],[0 1],2,2,2);
density_ref=zeros(2,2,2); density_ref(1,1,1)=1; density_ref(2,1,1)=1;
result=test_close(result,'xyz2pd interior bins',density,density_ref,1e-15,1e-15,...
                  'two in-range points must be counted in their x-axis bins and the out-of-range point discarded');

% Check nearest-neighbour lookup by Cartesian coordinate distance
[spin_idx,dist]=nearest_spin(spin_system,1);
result=test_close(result,'nearest_spin index',spin_idx,3,1e-15,1e-15,...
                  'spin three is the nearest neighbour of spin one in the fixture coordinates');
result=test_close(result,'nearest_spin distance',dist,0.5,1e-15,1e-15,...
                  'the nearest-neighbour distance is one half Angstrom');

% Check chemical-substance ownership lookup
result=test_close(result,'which_subst first part',which_subst(spin_system,[1 2]),1,1e-15,1e-15,...
                  'spins one and two both belong to the first chemical part');
result=test_close(result,'which_subst second part',which_subst(spin_system,3),2,1e-15,1e-15,...
                  'spin three belongs to the second chemical part');

% Check coupling extraction from forward and backward tensor entries
coupling=get_coupling(spin_system,1,2);
result=test_close(result,'get_coupling sum',coupling,diag([5 7 9]),1e-15,1e-15,...
                  'get_coupling must add the forward and backward coupling tensor cells');

% Check chemical shifts and Zeeman offsets from isotropic tensor parts
[cs_ppm,cs_hz]=chemshifts(spin_system);
result=test_close(result,'chemshifts ppm',[cs_ppm(1) cs_ppm(2) cs_ppm(3)],[1 2 0],1e-9,1e-12,...
                  'ppm shifts are one million times isotropic offsets divided by base frequencies');
result=test_close(result,'chemshifts hertz',[cs_hz(1) cs_hz(2) cs_hz(3)],[-100 -50 0],1e-7,1e-12,...
                  'hertz shifts are minus isotropic angular offsets divided by 2*pi');
result=test_close(result,'offsetof spin one',offsetof(spin_system,1),-100,1e-7,1e-12,...
                  'offsetof returns the same hertz offset for an individual spin');

% Check g-tensor extraction from the stored Zeeman scaling tensor
spin_system.inter.zeeman.ddscal{1}=2*eye(3);
g_ref=-spin_system.inter.zeeman.ddscal{1}*spin_system.inter.gammas(1)*...
      spin_system.tols.hbar/spin_system.tols.muB;
result=test_close(result,'gtensorof stored tensor',gtensorof(spin_system,1),g_ref,1e-15,1e-15,...
                  'gtensorof must apply the documented gamma*hbar/muB scaling to ddscal');

% Check replacement of isotropic tensor components without changing anisotropy
input_tensors={diag([1 2 6]),eye(3)};
shifted=shift_iso(input_tensors,1,10);
result=test_close(result,'shift_iso diagonal anisotropy',shifted{1},diag([8 9 13]),1e-12,1e-12,...
                  'replacing an isotropic part of three by ten adds seven to each diagonal component');
result=test_close(result,'shift_iso untouched tensor',shifted{2},eye(3),1e-15,1e-15,...
                  'tensors whose spin numbers are not requested must remain unchanged');

end


function spin_system=local_geometry_spin_system()

% Create quiet system output settings
spin_system.sys.output='hush';
spin_system.sys.disable={};

% Define spin identities and labels
spin_system.comp.nspins=3;
spin_system.comp.isotopes={'1H','13C','15N'};
spin_system.comp.labels={'proton','carbon','nitrogen'};
spin_system.comp.mults=[2 2 2];

% Define Cartesian coordinates and chemical parts
spin_system.inter.coordinates={[0 0 0],[2 0 0],[0.5 0 0]};
spin_system.chem.parts={[1 2],3};

% Define coupling tensor cells
spin_system.inter.coupling.matrix=cell(3,3);
spin_system.inter.coupling.matrix{1,2}=diag([1 2 3]);
spin_system.inter.coupling.matrix{2,1}=diag([4 5 6]);

% Define Zeeman frequencies with known isotropic offsets
basefrqs=2*pi*[100e6 25e6 10e6];
offsets=2*pi*[100 50 0];
spin_system.inter.basefrqs=basefrqs;
spin_system.inter.zeeman.matrix={basefrqs(1)*eye(3)+offsets(1)*eye(3),...
                                 basefrqs(2)*eye(3)+offsets(2)*eye(3),...
                                 basefrqs(3)*eye(3)+offsets(3)*eye(3)};

% Define magnetogyric ratios and fundamental constants
spin_system.inter.gammas=[spin('1H') spin('13C') spin('15N')];
spin_system.inter.zeeman.ddscal={eye(3),eye(3),eye(3)};
spin_system.tols.hbar=1.054571628e-34;
spin_system.tols.muB=9.274009994e-24;

end


