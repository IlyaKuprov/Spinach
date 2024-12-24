% Convergence and accuracy test of the spatial dynamics during a pulse
% field gradientspin (PFG) echo sequence. The accuracy spatial diffusi-
% on is tested with respect to the grid size and the finite difference
% stancil size. The PFG spin echo pulse sequence is simulated, and then
% the diffusion coefficient extracted back by fitting the Stejskal-Tan-
% ner equation.
%
% Run time: minutes on NVidia Tesla A100, much longer on CPU
%
% m.g.concilio@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function conv_test()

% Magnetic field
sys.magnet=11.7426; 

% Isotopes 
sys.isotopes={'1H'};

% Chemical shift
inter.zeeman.scalar={4.6};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.disable={'pt','krylov'};
sys.enable={'greedy','gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Aquisition parameters
parameters.sweep=5000; 
parameters.npoints=1024; 
parameters.spins={'1H'};
parameters.zerofill=32768;
parameters.axis_units='ppm';
parameters.offset=2500;      

% No relaxation
parameters.rlx_ph={};
parameters.rlx_op={};
          
% Gradient duration
parameters.delta_sml=0.002;

% Diffusion delay
parameters.delta_big=0.050;
  
% Sample length
parameters.dims=0.015;       % m  

% Diffusion coefficient
parameters.diff=18e-10;      % m^2/s

% Point stencil size list
stencil_sizes=[3 5 7];

% Gradient amplitude list 
grad_amps=linspace(0,0.5,20); % T/m

% Grid size list 
grid_sizes=ceil(linspace(1000,10000,50));

% Preallocate result arrays
D=zeros(numel(grid_sizes),numel(stencil_sizes));
T=zeros(numel(grid_sizes),numel(stencil_sizes));
sig_ints=cell(numel(grid_sizes),numel(stencil_sizes));

% Grid size loop
for k=1:numel(grid_sizes)
    
    % Set grid size 
    parameters.npts=grid_sizes(k);
    
    % Velocity field is zero
    parameters.u=zeros(parameters.npts,1);
    
    % Initial state - Gaussian in the middle
    parameters.rho0_ph={exp(-linspace(-10,10,parameters.npts).^2)'};
    parameters.rho0_st={state(spin_system,'Lz','1H','cheap')};
    
    % Detection state - uniform across sample
    parameters.coil_ph={ones(parameters.npts,1)};
    parameters.coil_st={state(spin_system,'L+','1H','cheap')};
    
    % Hush the reporting
    report(spin_system,'Spinach output hushed.');
    spin_system.sys.output='hush';
    
    % Stencil size loop
    for s=1:numel(stencil_sizes)
        
        % Set stencil size                                      
        parameters.deriv={'period',stencil_sizes(s)};  
        
        % Preallocate result array
        intensities=zeros(1,numel(grad_amps)); tic();
        
        % Gradient amplitude loop
        parfor n=1:numel(grad_amps)
            
            % Localise parameters
            localpar=parameters;
            
            % Set gradient amplitude
            localpar.g_amp=grad_amps(n);
            
            % Call Stejskal-Tanner sequence
            intensities(n)=imaging(spin_system,@st_ideal,localpar);
            
        end
        
        % Normalise and store the intensities
        sig_ints{k,s}=intensities/intensities(1);
        
        % Run S-T fitting (linear least squares via svd)
        delta_sml=parameters.delta_sml; delta_big=parameters.delta_big;
        expfactors=-(grad_amps*spin(parameters.spins{1})*delta_sml).^2*...
                    (delta_big-delta_sml/3)*1e-10;
        D(k,s)=log(sig_ints{k,s}(2:end))/expfactors(2:end); T(k,s)=toc();
        
        % Report to the user
        disp(['grid size ' num2str(grid_sizes(k)) ...
              ', stencil size ' num2str(stencil_sizes(s)) ...
              ', D=' num2str(D(k,s)) ', T=' num2str(T(k,s))]);
        
    end
    
end
    
% Generate convergence plot A
figure(); scale_figure([2.0 1.5]);
subplot(2,1,1); plot(grid_sizes',D); 
kgrid; ylim([0 20]); xlim([1000 10000]);
true_val=refline([0 18]);
set(true_val,'Color','k','LineStyle','--');
kxlabel('grid point count');
kylabel('$D / 10^{-10} m^2/s$');
klegend({'3-point stencil','5-point stencil',...
         '7-point stencil'},'Location','SouthEast');
subplot(2,1,2); plot(grid_sizes',T); 
kgrid; axis tight;
kxlabel('grid point count');
kylabel('wall clock time / s');
klegend({'3-point finite difference','5-point finite difference',...
         '7-point finite difference'},'Location','NorthWest');

% Generate convergence plot B
figure(); scale_figure([2.0 1.5]);
grids_to_plot=3:10;
for s=1:numel(stencil_sizes)
    subplot(2,3,s); hold on;
    for k=grids_to_plot
        plot(grad_amps,sig_ints{k,s});
    end
    box on; kgrid; xlim tight;
end
subplot(2,3,1); ktitle('3-point finite difference');
kylabel('normalised intensity');
subplot(2,3,2); ktitle('5-point finite difference');
subplot(2,3,3); ktitle('7-point finite difference');
ideal_curve=exp(-(grad_amps*spin(parameters.spins{1})*delta_sml).^2*...
                 (delta_big-delta_sml/3)*18e-10);
for s=1:numel(stencil_sizes)
    subplot(2,3,3+s); hold on;
    for k=grids_to_plot
        plot(grad_amps,abs(sig_ints{k,s}-ideal_curve));
        ylim([1e-3 1]); set(gca,'yscale','log');
    end
    box on; kgrid; xlim tight;
end
subplot(2,3,4); kylabel('difference from exact');
kxlabel('gradient amplitude, T/m');
subplot(2,3,5); kxlabel('gradient amplitude, T/m');
subplot(2,3,6); kxlabel('gradient amplitude, T/m');
subplot(2,3,3); klegend({[num2str(grid_sizes(3))  ' points'],...
                         [num2str(grid_sizes(4))  ' points'],...
                         [num2str(grid_sizes(5))  ' points'],...
                         [num2str(grid_sizes(6))  ' points'],...
                         [num2str(grid_sizes(7))  ' points'],...
                         [num2str(grid_sizes(8))  ' points'],...
                         [num2str(grid_sizes(9))  ' points'],...
                         [num2str(grid_sizes(10)) ' points']},...
                         'Location','SouthWest');

end

