% Numerical integral route to the Redfield relaxation superopera-
% tor. Syntax:
%
%            [R,dR]=ngce(spin_system,H0,H1,dt,tau_est,reg)
%
% Parameters:
%
%  H0 - static laboratory frame Hamiltonian commutation su-
%       peroperator acting in the background, a matrix
%
%  H1 - stochastic part (zero mean) of the laboratory frame
%       Hamiltonian commutation superoperator, a cell array
%       of matrices for each point in the MD trajectory.
%
%  dt - time step of the MD trajectory, seconds
%
%  tau_est - H1 autocorrelation time estimate for internal 
%            safety control, seconds
%
%  reg - optional overall relaxation rate, this is added to 
%        every eigenvalue of the resulting matrix to prevent
%        very small relaxation rates (e.g. singlets) from 
%        jumping into positive due to integration accuracy
%        limits and then causing problems
%           
% Outputs:
%
%  R  - laboratory frame relaxation superoperator
%
%  dR - standard deviation of the mean of R, element by element
%
% Note: enough trajectory points must be present to converge 
%       the ensemble averages and Redfield's integral.
%
% Note: the result is returned in the LABORATORY FRAME - eli-
%       minating non-secular terms is user's responsibility.
%
% ilya.kuprov@weizmann.ac.il
% jpresteg@uga.edu
%
% <https://spindynamics.org/wiki/index.php?title=ngce.m>

function [R,dR]=ngce(spin_system,H0,H1,dt,tau_est,reg)

% Check consistency
grumble(H0,H1,dt,tau_est);

% Coherent dynamics timescale
timescale=2*pi/normest(H0);

% Number of points in tau_c
npts_in_tau_c=ceil(tau_est/dt);

% Number of points under tau integral
n_tau_int_steps=5*npts_in_tau_c;

% Trajectory point count
traj_npts=numel(H1);

% Trajectory duration
traj_dur=dt*(traj_npts-1);

% Print timing diagnostics
report(spin_system,' ');
report(spin_system,['Shortest period in H0:           ' num2str(timescale,'%12.5e') ' seconds']);
report(spin_system,['Trajectory step length:          ' num2str(dt,'%12.5e') ' seconds']);
report(spin_system,['Total trajectory length:         ' num2str(traj_dur,'%12.5e') ' seconds']);
report(spin_system,['Total trajectory points:         ' num2str(traj_npts)]);
report(spin_system,['User estimate for tau_c:         ' num2str(tau_est,'%12.5e') ' seconds']);

% Enforce sufficient H0 sampling
report(spin_system,['Trajectory points per H0 period: ' num2str(timescale/dt)]);
if timescale/dt<50
    error('insufficient H0 dynamics sampling, reduce trajectory time step.');
end

% Enforce sufficient tau_c sampling
report(spin_system,['Trajectory points per tau_est:   ' num2str(tau_est/dt)]);
if tau_est/dt<10
    error('insufficient correlation function sampling, reduce trajectory time step.');
end

% Enforce sufficient trajectory duration
report(spin_system,['Trajectory duration / tau_est:   ' num2str(dt*(numel(H1)-1)/tau_est)]);
if traj_dur/tau_est<200
    error('insufficient ensemble average sampling, increase trajectory duration.');
end

% Get propagators for tau integrals
report(spin_system,'computing H0 exponentials...');
P=cell(n_tau_int_steps,1); P{1}=speye(size(H0));
P_dt=propagator(spin_system,H0,dt); 
for n=2:n_tau_int_steps % Keep cleaning up insignificant non-zeroes
    P{n}=clean_up(spin_system,P_dt*P{n-1},spin_system.tols.prop_chop);
end

% Determine the number of stripes
nstripes=floor(traj_npts/n_tau_int_steps);
report(spin_system,['trajectory cut into ' ...
       num2str(nstripes) ' ensemble instances']);

% Cut the trajectory up into stripes
H1=H1(1:nstripes*n_tau_int_steps);
H1=reshape(H1,[n_tau_int_steps nstripes]);

% Get unit state projector
U=unit_state(spin_system); UU=U*U';

% Compute trajectory stripe integrals
report(spin_system,'computing Redfield''s integral...');
calc_err=nargout>1;
rows=cell(nstripes,1); cols=cell(nstripes,1);
vals=cell(nstripes,1); vals_sq=cell(nstripes,1);
parfor s=1:nstripes % Inefficient, investigate
    
    % Get the stripe started
    H1s=H1(:,s); Ps=P; F=sparse(0);
    
    % Trapezium rule for tau integral
    f_curr=clean_up(spin_system,-dt*H1s{1}*Ps{1}*H1s{1}*Ps{1}',...
                    spin_system.tols.liouv_zero);
    for tau=2:n_tau_int_steps
        f_next=clean_up(spin_system,-dt*H1s{1}*Ps{tau}*H1s{tau}*Ps{tau}',...
                        spin_system.tols.liouv_zero);
        F=F+(f_curr+f_next)/2; f_curr=f_next;
    end

    % Keep real symmetric trace-preserving part
    F=real(F+F')/2; F=F-(U'*F*U)*UU;
    
    % Get stripe integral as a sparse array
    [stripe_rows,stripe_cols,stripe_vals]=find(F);
    rows{s}=stripe_rows; cols{s}=stripe_cols; vals{s}=stripe_vals;

    % Get stripe integral squared element by element
    if calc_err
        vals_sq{s}=stripe_vals.^2;
    end
    
end

% Assemble the superoperator and keep real symmetric part
rows=cell2mat(rows); cols=cell2mat(cols); vals=cell2mat(vals);
R_sum=sparse(rows,cols,vals,size(H0,1),size(H0,2));
R=R_sum/nstripes;

% Assemble standard deviation of the mean
if calc_err
    vals_sq=cell2mat(vals_sq);
    R_sq_sum=sparse(rows,cols,vals_sq,size(H0,1),size(H0,2));
    dR_var=(R_sq_sum-(R_sum.^2)/nstripes)/(nstripes*(nstripes-1));
    [rows,cols,vals]=find(dR_var);
    dR=sparse(rows,cols,sqrt(max(real(vals),0)),size(H0,1),size(H0,2));
end

% Apply regularisation
if exist('reg','var')&&(reg~=0)
    R=R-reg*unit_oper(spin_system);
end

% Make sure the unit state is not damped
R=R-(U'*R*U)*UU;

end

% Consistency enforcement
function grumble(H0,H1,dt,tau_est)
if ~isnumeric(H0)
    error('H0 must be a matrix.');
end
if ~iscell(H1)
    error('H1 must be a cell array of matrices.');
end
if (~isnumeric(dt))||(~isreal(dt))||(dt<=0)
    error('dt must be a positive real number.');
end
if (~isnumeric(tau_est))||(~isreal(tau_est))||(tau_est<=0)
    error('tau_est must be a positive real number.');
end
end

% "To engage in hand-to-hand combat, a trooper must have lost
% somewhere on the battlefield his rifle, pistol, knife, belt, 
% shovel, vest, and helmet. He must also have found a comple-
% tely clean spot with not a single stick or stone on it. And
% only then can he engage in ferocious hand-to-hand combat
% with a similar dumbass."
%
% A Russian hand-to-hand combat instructor

