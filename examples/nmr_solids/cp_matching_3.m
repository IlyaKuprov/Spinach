% Hartmann-Hahn matching condition test for a cross-polarisation
% experiment between a proton and a 15N nucleus. A 2D scan of
% power levels at a specific spinning rate.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function cp_matching_3()

% System specification
sys.magnet=9.394;
sys.isotopes={'1H','15N'};
          
% Interactions
inter.zeeman.scalar={0.1495  0.0000};
inter.coordinates={[-1.11551509    1.65289357   -1.19927242]
                   [-2.67552180    0.95825426    0.00000000]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.tols.inter_cutoff=5.0;
sys.tols.prox_cutoff=4.0;
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Relevant operators
Nx=operator(spin_system,'Lx','15N');
Ny=operator(spin_system,'Ly','15N');
Hx=operator(spin_system,'Lx','1H');
Hy=operator(spin_system,'Ly','1H');

% Power levels
powers=linspace(0e3,50e3,50);

% Experiment parameters
parameters.rate=10000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=3;
parameters.spins={'1H','15N'};
parameters.irr_opers={Hx Nx};
parameters.exc_opers={0*Hy 0*Ny};
parameters.rho0=state(spin_system,'Lx','1H');
parameters.coil=state(spin_system,'Lx','15N');
parameters.time_steps=4e-5*ones(1,10);
parameters.grid='rep_2ang_200pts_oct';

% Parallel loop over power levels
cp=zeros(numel(powers),numel(powers)); figure();
for n=1:numel(powers)
    
    parfor k=1:numel(powers)
        
        % MAS parameters
        localpar=parameters;
        localpar.irr_powers=[powers(n)*ones(1,10); 
                             powers(k)*ones(1,10)]; %#ok<PFBNS>
        
        % Simulation
        fid=singlerot(spin_system,@cp_contact_hard,localpar,'nmr');
        cp(n,k)=real(fid(end));
        
    end
    
    % Plot as an image
    imagesc([min(powers) max(powers)],...
            [min(powers) max(powers)],cp);
    kxlabel('$^{1}$H spin-lock RF power, Hz'); 
    kylabel('$^{15}$N spin-lock RF power, Hz'); 
    set(gca,'YDir','normal'); drawnow();
    
end

end

