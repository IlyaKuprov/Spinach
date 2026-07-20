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

% Main magnet
sys.magnet=0.13;

% Eliminate dipolar hyperfine (will
% get recomputed from coordinates)
inter.coupling.matrix{1,2}=[];

% Create three such systems
[sys,inter]=merge_inp({sys,sys,sys},{inter,inter,inter});

% Specify Cartesian coordinates
inter.coordinates={[0  0 -50]; [0  0 -51.5];
                   [0  10  0]; [0  10  1.5];
                   [50 0   0]; [50 0   1.5]};

% Relax tolerances
sys.tols.prox_cutoff=Inf;   % Angstrom
sys.tols.inter_cutoff=1.0;  % Hz

% Modest parallel pool
sys.parallel={'processes',8};

% Spinach housekeeping
spin_system=create(sys,inter);
 
% Basis set and formalism
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=1;

% Spinach housekeeping
spin_system=basis(spin_system,bas);

% Set assumptions
spin_system=assume(spin_system,'esr');

% Make a Hamiltonian
[I,Q]=hamiltonian(spin_system);

% Specify an orientation
H=I+orientation(Q,[0 0 0]);

% Load coil data
load('grad_coils_data.mat','B_X_Coils',...
                           'B_Y_Coils',... 
                           'B_Z_Coils');

% Remove average Z field
B_Z_Coils=B_Z_Coils-mean(B_Z_Coils);

% Get magnetic fields at spin locations
B0G=coil_map_lookup([0 0 0],B_X_Coils,B_Y_Coils,B_Z_Coils,...
                    [0 0 100],inter.coordinates);

% Add Zeeman terms manually
for n=1:numel(sys.isotopes)

    % Get the operators
    Lx=operator(spin_system,{'Lx'},{n});
    Ly=operator(spin_system,{'Ly'},{n});
    Lz=operator(spin_system,{'Lz'},{n});

    % Get the frequencies
    omega_x=-spin_system.inter.gammas(n)*B0G{n}(1)*0;
    omega_y=-spin_system.inter.gammas(n)*B0G{n}(2)*0;
    omega_z=-spin_system.inter.gammas(n)*B0G{n}(3);

    % Add to the Hamiltonian
    H=H+omega_x*Lx+omega_y*Ly+omega_z*Lz;

end

% Set the sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E');
parameters.coil=state(spin_system,'L+','E');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=2e9;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.axis_units='GHz-labframe';

% Simulation
fid=evolution(spin_system,H,parameters.coil,parameters.rho0,...
              1/parameters.sweep,parameters.npoints-1,'observable');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Perform Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plot the spectrum
kfigure(); plot_1d(spin_system,real(spectrum),parameters);    

end

