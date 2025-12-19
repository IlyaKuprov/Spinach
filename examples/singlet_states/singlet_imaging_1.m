% Singlet imaging in a system with one-dimensional 
% diffusion and flow.
%
% Calculation time: minutes
%
% a.j.allami@soton.ac.uk
% g.pileio@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function singlet_imaging_1()

% Spin system and interactions
sys.magnet=9.4;
sys.isotopes={'13C','13C'};
inter.zeeman.scalar={0.03,-0.03};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=55;
inter.coordinates={[0.00 0.00 0.00];
                   [1.20 0.00 0.00]};

% Relaxation theory
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={1e-9};

% Relaxation superoperator accuracy
sys.tols.rlx_integration=1e-5;
sys.tols.rlx_zero=1e-5;

% Algorithmic options
sys.enable={'greedy'};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sample geometry
parameters.dims=[0.10 0.015];
parameters.npts=[150 15];
parameters.deriv={'period',7};

% Sequence parameters
parameters.spins={'13C'};
parameters.offset=0;

% Assumptions
spin_system=assume(spin_system,'nmr');

% Relaxation phantom
parameters.rlx_ph={ones(parameters.npts)};
parameters.rlx_op={relaxation(spin_system)};

% NMR sample tube phantom
tube=zeros(parameters.npts); tube(:,6:10)=1;

% Initial and detection state phantoms
parameters.rho0_ph={tube};
parameters.rho0_st={state(spin_system,'Lz','13C','cheap')};
parameters.coil_ph={tube};
parameters.coil_st={state(spin_system,'L+','13C','cheap')};

% Diffusion and flow
parameters.u=-6e-2*ones(parameters.npts);
parameters.v=zeros(parameters.npts);
parameters.diff=[3.6e-6 0 0; 0 0 0; 0 0 0];

    % Create a pulse sequence in the imaging context
    function traj=tube_flow(spin_system,parameters,H,R,K,G,F)
        
        % Compose Liouvillian
        L=H+F+1i*R+1i*K;
        
        % Get pulse operators
        Lx=operator(spin_system,'Lx',parameters.spins{1});
        Ly=operator(spin_system,'Ly',parameters.spins{1});
        Lx=kron(speye(prod(parameters.npts)),Lx);
        Ly=kron(speye(prod(parameters.npts)),Ly);

        % Pulse phase
        rf_phi=pi/2;

        % Number of steps in the pulse
        pulse_nsteps=25;

        % Overall pulse duration
        pulse_time=2e-3;

        % Pulse frequency
        pulse_frq=2000;

        % Pulse amplitude table
        rf_amp_list=2*pi*100*pulse_shape('gaussian',pulse_nsteps);

        % Pulse block duration table
        rf_dur_list=(pulse_time/pulse_nsteps)*ones(1,pulse_nsteps);

        % Gradient amplitude
        grad_amp=6e-3; % T/m

        % Apply a pulse at a specific frequency under Gx
        rho=shaped_pulse_af(spin_system,L+grad_amp*G{1},Lx,Ly,parameters.rho0,...
                            pulse_frq*ones(1,pulse_nsteps),rf_amp_list,rf_dur_list,...
                            rf_phi,2,'expv');

        % Keep magnetisation
        rho_magn=rho;
        
        % M2S parameters
        J=55; delta_v=6;
        t=1/(4*sqrt(J^2+delta_v^2));
        m1=floor(pi*J/(2*delta_v));
        if mod(m1,2)~=0, m1=m1+1; end
        
        % M2S sequence
        for k=1:m1
            rho=step(spin_system,H,rho,t);
            rho=step(spin_system,Lx,rho,pi);
            rho=step(spin_system,H,rho,t);
        end
        rho=step(spin_system,Lx,rho,pi/2);
        rho=step(spin_system,H,rho,t);
        for k=1:m1/2
            rho=step(spin_system,H,rho,t);
            rho=step(spin_system,Lx,rho,pi);
            rho=step(spin_system,H,rho,t);
        end
        rho_sing=-step(spin_system,H,rho,t);
        
        % Free evolution
        traj{1}=evolution(spin_system,L,[],rho_sing,0.01,100,'trajectory');
        traj{2}=evolution(spin_system,L,[],rho_magn,0.01,100,'trajectory');
        
    end

% Run the pulse sequence in the imaging context
traj=imaging(spin_system,@tube_flow,parameters);

% Get detection states
coil_sing=singlet(spin_system,1,2);
coil_magn=state(spin_system,'L+','all');
                
% Show the movie
kfigure();
for n=1:100
    tube_sing=fpl2phan(traj{1}(:,n),coil_sing,parameters.npts);
    tube_magn=fpl2phan(traj{2}(:,n),coil_magn,parameters.npts);
    subplot(1,2,1); imagesc(real(tube_sing),25*[-2e-3 6e-3]/2); 
    kcolourbar; ktitle('singlet state'); kxlabel('mm');
    subplot(1,2,2); imagesc(real(tube_magn),25*[-2e-3 6e-3]); 
    kcolourbar; ktitle('transv. magn.'); kxlabel('mm');
    drawnow; pause(0.01);
end

end

