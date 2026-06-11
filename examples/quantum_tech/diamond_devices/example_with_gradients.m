% Three-P1 Spinach Hamiltonian with measured gradient coil fields.
%
% Syntax:
%
%   [H,spin_system,gradient]=example_with_gradients()
%   [H,spin_system,gradient]=example_with_gradients(currents)
%
% Inputs:
%
%   currents - X/Y/Z gradient coil currents, A
%
% Outputs:
%
%   H           - Spinach Hamiltonian, rad/s
%   spin_system - Spinach spin system object
%   gradient    - gradient map metadata and interpolated fields

function [H,spin_system,gradient]=example_with_gradients()

% Get a single defect spin system
parameters.orientation='111';
parameters.nitrogen='15N';
[sys,inter]=diamond_p1(parameters);

% Eliminate dipolar hyperfine (will
% get recomputed from coordinates)
inter.coupling.matrix{1,2}=[];

% Create three such systems
[sys,inter]=merge_inp({sys,sys,sys},{inter,inter,inter});

% Specify Cartesian coordinates
inter.coordinates={[0 0 -50]; [0 0 -51.5];
                   [0 10  0]; [0 10 1.5];
                   [50 0  0]; [50 0 1.5]};

% Load coil data
load('grad_coils_data.mat','B_X_Coils',...
                           'B_Y_Coils',... 
                           'B_Z_Coils');

% Get magnetic fields at spin locations
B1=coil_map_lookup([0 0 0],B_X_Coils,B_Y_Coils,B_Z_Coils,...
                   [1.0 1.0 1.0],inter.coordinates);











% spin_op=operator(spin_system,'Lz',1);
% H_grad=sparse(size(spin_op,1),size(spin_op,2));
% 
% for n=1:spin_system.comp.nspins
%     zeeman_per_tesla=spin_system.inter.zeeman.matrix{n}/ ...
%         spin_system.inter.magnet;
%     omega=zeeman_per_tesla*fields(n,:).';
%     H_grad=H_grad+omega(1)*operator(spin_system,'Lx',n)+ ...
%         omega(2)*operator(spin_system,'Ly',n)+ ...
%         omega(3)*operator(spin_system,'Lz',n);
% end
% 
% 
% 
% 
% 
% 
% 
% fields=fields+currents(n)*interpolate_field_map(map,spin_coords_m);
% 
% 
% field=zeros(size(positions));
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % Load gradient coil fields
% load('grad_coils_data.mat',
% 
% gradient.data_file=fullfile(repo_root,'data','field_maps', ...
%     'surface_micro_gradients','grad_coils_data_21_5_26.mat');
% gradient.currents=[1 1 1]; % A, X/Y/Z gradient coil currents
% gradient.fields=gradient_fields_at_spins(gradient.data_file, ...
%     inter.coordinates,gradient.currents);
% 
% % Spinach parallel setup
% sys.parallel={'local',1};
% 
% % Cutoff modifications
% sys.tols.prox_cutoff=Inf;   % Angstrom
% sys.tols.inter_cutoff=1.0;  % Hz
% 
% 
% 
% % Spinach housekeeping
% spin_system=create(sys,inter);
% 
% % Basis set and formalism
% bas.formalism='zeeman-hilb';
% bas.approximation='none';
% 
% % Spinach housekeeping
% spin_system=basis(spin_system,bas);
% 
% % Set assumptions
% spin_system=assume(spin_system,'esr');
% 
% % Make a Hamiltonian
% [I,Q]=hamiltonian(spin_system);
% 
% % Specify an orientation
% H=I+orientation(Q,[0 0 0]);
% 
% % Add position-dependent Zeeman terms from gradient coil fields
% H=H+gradient_zeeman_hamiltonian(spin_system,gradient.fields);
% 
% end



end




