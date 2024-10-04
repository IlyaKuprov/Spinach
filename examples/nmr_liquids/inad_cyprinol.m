% INADEQUATE spectrum of cyprinol. The sequence filters out double-
% quantum coherences and thus only keeps pairs of spins that have
% a J-coupling. A parallel sum over isotopomers that have adjacent
% 13C spins is used.
%
% Calculation time: minutes
% 
% Bud Macaulay, Ilya Kuprov

function inad_cyprinol()

% Spin system - cyprinol
[sys,inter]=cyprinol();

% Magnet field
sys.magnet=11.7;

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.space_level=1;
bas.level=4;

% seq parameters
parameters.spins={'13C'};
parameters.J=50;
parameters.decouple={'1H'};
parameters.offset=5000;
parameters.sweep=10000;
parameters.npoints=4096;
parameters.zerofill=8192;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Spinach housekeeping
spin_system=create(sys,inter);

% Generate isotopomers
subsystems=dilute(spin_system,'13C',2);

% Preallocate the answer
spectrum=zeros([parameters.zerofill(1) 1],'like',1i);

% Loop over isotopomers
parfor n=1:numel(subsystems)
    
    % Check the J-coupling between the two carbons
    c_idx=find(cellfun(@(x)strcmp(x,'13C'),subsystems{n}.comp.isotopes)); 
    J=trace(get_coupling(subsystems{n},c_idx(1),c_idx(2)))/3;
    
    % For couplings stronger than 1 Hz
    if abs(J)>2*pi*1.0      
        
        % Build the basis
        subsystem=basis(subsystems{n},bas);
        
        % Simulation
        fid=liquid(subsystem,@inadequate,parameters,'nmr');
        
        % Apodisation
        fid=apodisation(spin_system,fid,{{'exp',6}});
        
        % Fourier transform
        spectrum=spectrum+fftshift(fft(fid,parameters.zerofill));
        
    end
    
end

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

