% Tests graph, geometry, lattice, and coordinate utilities. Syntax:
%
%                    result=test_graph_geometry_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks small graph decompositions, lattice generation,
% dihedral angles, coordinate density binning, nearest-spin lookup,
% substance and label lookup helpers, and coordinate-derived dipolar
% tensor identities.
%
% ilya.kuprov@weizmann.ac.il

function result=test_graph_geometry_suite()

% Announce the test target
fprintf('TESTING: Graph, geometry, and coordinate utilities\n');

% State the utility target of the test
result=new_test_result('kernel/graph_geometry_suite',...
                       'Graph, geometry, and coordinate utilities',...
                       'Small graph and geometry helpers must reproduce exact combinatorial and tensor references.');

% Check a two-period cubic lattice coordinate and periodic-cell layout
[sys,inter]=cubic_lattice('13C',1.5,2);
coord_obs=sortrows(cell2mat(inter.coordinates));
coord_ref=sortrows(1.5*[0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1]);
result=test_true(result,'cubic_lattice isotope count',numel(sys.isotopes)==8,...
                 'a two-period cubic lattice contains 2^3 isotope entries');
result=test_close(result,'cubic_lattice coordinates',coord_obs,coord_ref,...
                  1e-15,1e-15,...
                  'coordinates enumerate all corners of the requested cubic grid');
result=test_close(result,'cubic_lattice first pbc vector',inter.pbc{1},[3 0 0],...
                  1e-15,1e-15,...
                  'the first periodic translation is spacing*n_periods along x');

% Check a right-handed coordinate set with a ninety-degree dihedral
phi=dihedral([1 0 0],[0 0 0],[0 1 0],[0 1 1]);
result=test_close(result,'dihedral right angle magnitude',abs(phi),90,...
                  1e-13,1e-13,...
                  'orthogonal adjacent planes give a ninety-degree dihedral magnitude');

% Check depth-first partitioning on a three-node path graph
path_graph=logical([0 1 0;1 0 1;0 1 0]);
subgraphs=dfpt(path_graph,2);
subgraph_ref=logical([1 1 0;0 1 1]);
result=test_true(result,'dfpt path graph pairs',isequal(full(subgraphs),subgraph_ref),...
                 'size-two connected subgraphs of a three-node path are its two edges');

% Check strongly connected components on two disconnected two-cycles
strong_graph=logical([0 1 0 0;1 0 0 0;0 0 0 1;0 0 1 0]);
component_idx=scomponents(strong_graph);
component_ok=(component_idx(1)==component_idx(2))&&...
             (component_idx(3)==component_idx(4))&&...
             (component_idx(1)~=component_idx(3));
result=test_true(result,'scomponents two cycles',component_ok,...
                 'two disconnected two-cycles form two distinct strongly connected components');

% Check three-dimensional point-density binning with one discarded point
points=[0.25 0.25 0.25;0.75 0.75 0.75;1.25 0.25 0.25];
density_obs=xyz2pd(points,[0 1],[0 1],[0 1],2,2,2);
density_ref=zeros(2,2,2);
density_ref(1,1,1)=1;
density_ref(2,2,2)=1;
result=test_close(result,'xyz2pd grid counts',density_obs,density_ref,...
                  1e-15,1e-15,...
                  'points inside the grid are counted in their cells and outside points are discarded');

% Check nearest-spin lookup from Cartesian coordinates
spin_system=local_geometry_system();
[near_idx,near_dist]=nearest_spin(spin_system,1);
result=test_close(result,'nearest_spin index',near_idx,3,1e-15,1e-15,...
                  'the closest coordinate to spin one is spin three');
result=test_close(result,'nearest_spin distance',near_dist,0.5,1e-15,1e-15,...
                  'the nearest-spin distance is the Euclidean coordinate separation');

% Check label and chemical-substance lookup helpers
sys.labels={'ha','ca','hb'};
sys.isotopes={'1H','13C','1H'};
result=test_close(result,'idxof label lookup',idxof(sys,'ca'),2,1e-15,1e-15,...
                  'a unique spin label resolves to its isotope-list index');
result=test_close(result,'which_subst same substance',which_subst(spin_system,[1 3]),1,...
                  1e-15,1e-15,...
                  'spins one and three belong to the first chemical substance');

% Check coordinate-derived dipolar coupling tensor identities on the z axis
[dip_coupling,alpha,beta,gamma,dip_mat]=xyz2dd([0 0 0],[0 0 1],'1H','1H');
result=test_close(result,'xyz2dd z-axis Euler angles',[alpha beta gamma],[0 0 0],...
                  1e-13,1e-13,...
                  'a z-axis dipolar vector has zero alpha, beta, and gamma in this convention');
result=test_close(result,'xyz2dd axial tensor',dip_mat,dip_coupling*diag([1 1 -2]),...
                  1e-10,1e-12,...
                  'a z-axis point dipole gives the axial traceless dipolar tensor');
[dip_half,~,~,~]=xyz2dd([0 0 0],[0 0 2],'1H','1H');
result=test_close(result,'xyz2dd inverse-cube scaling',dip_half,dip_coupling/8,...
                  1e-8,1e-12,...
                  'doubling the internuclear distance divides the dipolar coupling by eight');

% Check coordinate-derived hyperfine tensor symmetry and electron-count scaling
hfc_one=xyz2hfc([0 0 0],[0 0 1],'1H',1);
hfc_two=xyz2hfc([0 0 0],[0 0 1],'1H',2);
result=test_close(result,'xyz2hfc symmetry',hfc_one,hfc_one',1e-12,1e-12,...
                  'the point-dipole hyperfine tensor is symmetric');
result=test_close(result,'xyz2hfc tracelessness',trace(hfc_one),0,1e-12,1e-12,...
                  'the point-dipole hyperfine tensor is traceless');
result=test_close(result,'xyz2hfc electron-count scaling',hfc_two,2*hfc_one,...
                  1e-12,1e-12,...
                  'the point-dipole hyperfine tensor scales linearly with the number of electrons');

end

% Minimal spin_system for geometry and chemical-substance helpers
function spin_system=local_geometry_system()

spin_system.comp.nspins=3;
spin_system.comp.isotopes={'1H','13C','1H'};
spin_system.chem.parts={[1 3],2};
spin_system.inter.coordinates={[0 0 0],[2 0 0],[0.5 0 0]};

end


