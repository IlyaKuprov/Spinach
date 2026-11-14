% Deterministic lattice, occupancy, geometry and frame tests.
%
%                         result=test_cce_bath()
%
% ilya.kuprov@weizmann.ac.il

function result=test_cce_bath()

% Initialize standard harness fields only
result=new_test_result('kernel/cce_bath','CCE bath construction',...
                       'Nonorthogonal cells, occupancy, density and covariance.');

% Compare nonorthogonal replication against direct Cartesian arithmetic
s=struct('mode','lattice','isotopes',{{'1H','13C'}},'magnetic',1,...
         'cell_A',[2 0 0;1 3 0;0.2 0.4 4],'fractional',[0 0 0;0.25 0.5 0.75],...
         'repeats',[2 1 1]);
[sys,~,bath]=cce_bath(s); expected=[s.fractional*s.cell_A;(s.fractional+[1 0 0])*s.cell_A];
result=test_close(result,'Nonorthogonal row-vector cell convention',bath.coordinates_A,expected,0,0,'Nonorthogonal row-vector cell convention');
result=test_true(result,'Unique physical site IDs',numel(sys.isotopes)==4&&numel(unique(bath.ids))==4,'Unique physical site IDs');

% Check isotope occupancy without mutating the global random stream
s.species={'1H','2H'}; s.occupancy=[0.4 0.6]; s.seed=791;
state=rng; [~,~,a]=cce_bath(s); [~,~,b]=cce_bath(s);
result=test_true(result,'Seeded isotope occupancy is reproducible and local',isequal(a,b)&&isequal(state,rng),'Seeded isotope occupancy is reproducible and local');

% Verify active frame rotation of coordinates and Cartesian interactions
q=struct('mode','explicit','isotopes',{{'1H','13C'}},'magnetic',1,...
         'coordinates_A',[1 2 3;4 5 6],'rotation',[0 -1 0;1 0 0;0 0 1],...
         'translation_A',[2 3 4]);
q.pairs_Hz=struct('ids',[1 2],'tensor',[1 2 3;4 5 6;7 8 10]);
[~,inter,bath]=cce_bath(q);
result=test_close(result,'Coordinate covariance',bath.coordinates_A,q.coordinates_A*q.rotation.'+q.translation_A,0,0,'Coordinate covariance');
result=test_close(result,'Tensor covariance',inter.coupling.matrix{1,2},q.rotation*q.pairs_Hz.tensor*q.rotation.',0,0,'Tensor covariance');

% Test a two-dimensional number density and finite packing exclusion
r=struct('mode','random','isotopes',{{'E'}},'magnetic',0,...
         'bounds_A',[0 0 0;100 200 0],'dimension',2,'density',1e17,...
         'seed',19,'max_attempts',10000,'exclusion_A',1);
[~,~,bath]=cce_bath(r);
result=test_true(result,'Surface density gives twenty unique sites',numel(bath.ids)==20&&size(unique(bath.coordinates_A,'rows'),1)==20,'Surface density gives twenty unique sites');
result=test_close(result,'SI density and area conversion',bath.realized_density,r.density,1,1e-14,'SI density and area conversion');

% Reject duplicate positions and remove sites before preserving group closure
q=rmfield(q,{'rotation','translation_A','pairs_Hz'}); q.coordinates_A=[0 0 0;1 0 0];
q.groups={[1 2]}; q.remove_ids=2; [~,~,bath]=cce_bath(q);
result=test_true(result,'Defect removal preserves physical group IDs',isequal(bath.ids,1)&&isequal(bath.groups,{1}),'Defect removal preserves physical group IDs');
q=rmfield(q,'remove_ids'); q.coordinates_A=zeros(2,3); rejected=false;
try, cce_bath(q); catch, rejected=true; end
result=test_true(result,'Duplicate physical sites are refused',rejected,'Duplicate physical sites are refused');

% Check explicit physical grouping and isotope row normalization
q=struct('mode','explicit','isotopes',{{'1H';'13C'}},'magnetic',1,...
         'coordinates_A',[0 0 0;1 0 0],'groups',{{1,2}});
[sys,inter,b]=cce_bath(q);
result=test_true(result,'Complete groups and isotope row convention',...
                 isequal(b.groups,{1,2})&&isrow(sys.isotopes),...
                 'User groups explicitly partition every retained site.');

% Check scalar isotropic interactions match explicitly supplied tensors
q.pairs_Hz=struct('ids',[1 2],'tensor',37); [~,scalar]=cce_bath(q);
q.pairs_Hz.tensor=37*eye(3); [~,tensor]=cce_bath(q);
result=test_close(result,'Scalar isotropic tensor',scalar.coupling.matrix{1,2},...
                  tensor.coupling.matrix{1,2},0,0,'Both inputs are in Hz.');

% Check bath-to-Spinach point-dipole covariance without a Hamiltonian shortcut
q=rmfield(q,'pairs_Hz'); R=[0 -1 0;1 0 0;0 0 1];
q.rotation=R; q.translation_A=[3 4 5];
[rsys,rinter]=cce_bath(q); sys.output='hush'; rsys.output='hush';
sys.parallel={'local',1}; rsys.parallel={'local',1};
plain=create(sys,inter); rotated=create(rsys,rinter);
T=plain.inter.coupling.matrix{1,2}; Tr=rotated.inter.coupling.matrix{1,2};
result=test_close(result,'Native dipolar transform covariance',Tr,R*T*R.',...
                  100*eps*norm(T,'fro'),100*eps,'Rigid translations leave dipoles invariant.');

% Reject invalid original site labels rather than silently dropping them
q=rmfield(q,{'rotation','translation_A'}); q.groups={[1 3]};
rejected=false; try, cce_bath(q); catch, rejected=true; end
result=test_true(result,'Unknown group site refused',rejected,'Original candidate IDs are checked.');
q.groups={1,2}; q.pairs_Hz=struct('ids',[1 3],'tensor',eye(3));
rejected=false; try, cce_bath(q); catch, rejected=true; end
result=test_true(result,'Unknown pair site refused',rejected,'A misspelled interaction cannot disappear.');

% Check the linear density boundary independently of the areal case
r=struct('mode','random','isotopes',{{'1H'}},'magnetic',0,'seed',12,...
         'dimension',1,'bounds_A',[0 2 3;10 2 3],'density',3e9,'max_attempts',20);
[~,~,line]=cce_bath(r);
result=test_true(result,'Linear density and frozen coordinates',numel(line.ids)==3&&...
                 all(line.coordinates_A(:,2)==2)&&all(line.coordinates_A(:,3)==3),...
                 'Three spins in 10 Angstrom at 3e9 per metre.');

% Expose unattainable random packing at the caller ceiling
r.count=2; r=rmfield(r,'density'); r.exclusion_A=20; r.max_attempts=2;
rejected=false; try, cce_bath(r); catch e, rejected=strcmp(e.identifier,'cce_bath:packing'); end
result=test_true(result,'Packing ceiling exposed',rejected,'No hidden relaxation of exclusion radius.');

fprintf('CCE bath: lattice error %.3g; density %.16g m^-2; all geometric checks complete.\n',max(abs(expected-a.candidate_coordinates_A),[],'all'),2e18);

end

