% Cross-polarisation experiment in the rotating frame. Applies an
% ideal pi/2 pulse using the specified operators, then evolves the
% system with the specified spin-lock terms added to the Liovilli-
% an. The contact curve is returned. Syntax:
%
%   contact_curve=cp_contact_hard(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.irr_powers - a matrix containing the values
%                             of the spin-lock nutation fre-
%                             quency on each channel (rows)
%                             at each time slice (cols), Hz
%                            
%     parameters.irr_opers  - a cell array of spin operators
%                             corresponding to the spin-lock
%                             on each channel
%
%     parameters.exc_opers  - a cell array of spin operators
%                             for the ideal pi/2 excitation 
%                             pulse (same flip angle on all
%                             channels)
%
%     parameters.npoints    - Number of Ticks the Rotor Makes                             
%
%     parameters.rho0       - initial state vector
%
%     parameters.coil       - detection state vector
%
%     H - Hamiltonian matrix, received from context function
%
%     R - relaxation superoperator, received from context function
%
%     K - kinetics superoperator, received from context function
%
% Output:
%
%     contact_curve - contact curve detected on the coil
%                     state specified in parameters.coil
%
% marina.carravetta@soton.ac.uk
% p.t.williamson@soton.ac.uk
% guinevere.mathies@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
% m.keitel@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=cp_contact_hard.m>

function contact_curve=cp_contact_hard_tdhamil(spin_system,parameters,H,R,K)

% Get pulse sequence timing
rotor_period=1/parameters.rate;
nintspp=2*parameters.max_rank+1;
nsteps=nintspp*parameters.nperiods;
dt=rotor_period/nintspp;

% Preallocate the observable array
contact_curve=zeros([1 nsteps+1],'like',1i);

% Get the operators started
exc_oper=sparse(0); irr_oper=sparse(0);

% Build the excitation operator
for n=1:numel(parameters.exc_opers)
    exc_oper=exc_oper+parameters.exc_opers{n};
end

% Build the spin-lock operator
for k=1:numel(parameters.irr_opers)
    irr_oper=irr_oper+2*pi*parameters.irr_powers(k)*parameters.irr_opers{k};
end

% Apply a perfect pi/2 excitation pulse
rho=step(spin_system,exc_oper,parameters.rho0,pi/2);

% First point of the trajectory
contact_curve(:,1)=hdot(parameters.coil,rho);

% Loop over the time steps
for n=1:nsteps

    % Get the current evolution generator
    G=H{mod(n-1,nintspp)+1}+1i*R+1i*K+irr_oper;
    
    % Take the evolution step
    rho=step(spin_system,G,rho,dt);

    % Get the observable quantity
    contact_curve(:,n+1)=hdot(parameters.coil,rho);

end

end

