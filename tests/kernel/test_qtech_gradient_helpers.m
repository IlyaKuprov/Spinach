% Tests quantum technology gradient field helpers. Syntax:
%
%                    result=test_qtech_gradient_helpers()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% a.arnab@weizmann.ac.il

function result=test_qtech_gradient_helpers()

% Announce the test target
fprintf('TESTING: quantum technology gradient helpers\n');

% State the target of the test
result=new_test_result('kernel/qtech_gradient_helpers',...
                       'Quantum technology gradient helpers',...
                       'field-map interpolation and local Zeeman Hamiltonian assembly must be consistent.');

% Build simple affine field maps on a cube
[x,y,z]=ndgrid([0 1],[0 1],[0 1]);
map_a=[x(:) y(:) z(:) ...
       2*x(:)+3*y(:)+5*z(:) ...
      -1*x(:)+7*y(:) ...
       11*z(:)];
map_b=[x(:) y(:) z(:) ...
      -3*x(:)+2*y(:)+z(:) ...
       5*x(:)-y(:)+4*z(:) ...
       7*x(:)];
field_maps.coil_a=map_a;
field_maps.coil_b=map_b;

% Interpolate them between grid points
positions=[0.25 0.50 0.75; 0.75 0.25 0.50];
fields=gradient_fields_at_spins(field_maps,positions,[2 -1],{'coil_a','coil_b'});
fields_ref=2*[2*positions(:,1)+3*positions(:,2)+5*positions(:,3) ...
             -positions(:,1)+7*positions(:,2) ...
              11*positions(:,3)]-...
           [-3*positions(:,1)+2*positions(:,2)+positions(:,3) ...
             5*positions(:,1)-positions(:,2)+4*positions(:,3) ...
             7*positions(:,1)];
result=test_close(result,'field map interpolation',fields,fields_ref,...
                  1e-12,1e-12,'linear interpolation must reproduce affine maps exactly.');

% Check Spinach coordinate cell input
coord_cells={positions(1,:)*1e10,positions(2,:)*1e10};
fields=gradient_fields_at_spins(field_maps,coord_cells,[2 -1],{'coil_a','coil_b'});
result=test_close(result,'coordinate cell conversion',fields,fields_ref,...
                  1e-12,1e-12,'Spinach coordinate cells must be converted from Angstrom to metres.');

% Check field map loading from MAT files
map_file=[tempname '.mat'];
file_cleaner=onCleanup(@()delete(map_file));
save(map_file,'-struct','field_maps');
fields=gradient_fields_at_spins(map_file,positions,[2 -1],{'coil_a','coil_b'});
result=test_close(result,'MAT file field maps',fields,fields_ref,...
                  1e-12,1e-12,'field maps must be accepted from MAT files by name.');

% Build a one-proton Hilbert-space spin system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={1};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Compare helper output to explicit Cartesian Zeeman construction
local_field=[1e-6 2e-6 3e-6];
H_obs=gradient_zeeman_hamiltonian(spin_system,local_field);
zeeman_per_tesla=spin_system.inter.zeeman.matrix{1}/spin_system.inter.magnet;
omega=zeeman_per_tesla*local_field.';
H_ref=omega(1)*operator(spin_system,'Lx',1)+...
      omega(2)*operator(spin_system,'Ly',1)+...
      omega(3)*operator(spin_system,'Lz',1);
result=test_close(result,'local Zeeman Hamiltonian',H_obs,H_ref,...
                  1e-12,1e-12,'local fields must multiply Cartesian spin operators.');

% Compare axial field output to the existing uniform Zeeman pathway
local_field=[0 0 3e-6];
H_obs=gradient_zeeman_hamiltonian(spin_system,local_field);
H_ref=local_field(3)*hamiltonian(assume(spin_system,'labframe','zeeman'))/...
      spin_system.inter.magnet;
result=test_close(result,'uniform axial Zeeman path',H_obs,H_ref,...
                  1e-12,1e-12,'uniform axial fields must match the existing Zeeman Hamiltonian pathway.');

end

