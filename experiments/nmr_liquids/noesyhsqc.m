% Phase-sensitive NOESY-HSQC sequence from
%
%         http://pubs.acs.org/doi/abs/10.1021/bi00441a004 
%
% using the bidirectional propagation method described in
%
%          http://dx.doi.org/10.1016/j.jmr.2014.04.002
%
% The sequence is hard-wired to {F1,F2,F3}={1H,15N,1H} with carbon
% decoupled throughout. Syntax:
%
%          fid=noesyhsqc(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.npoints   - a vector of three integers giving the
%                           number of points in the three temporal
%                           dimensions, ordered as [t1 t2 t3].
%
%    parameters.sweep     - a vector of three real numbers giving
%                           the sweep widths in the three frequen-
%                           cy dimensions, ordered as [f1 f2 f3].
%
%    parameters.J         - J-coupling for the HSQC stage magneti-
%                           sation transfer, Hz.
%
%    parameters.tmix      - NOESY stage mixing time, seconds.
%
%    H - Hamiltonian matrix, received from context function
%
%    R - relaxation superoperator, received from context function
%
%    K - kinetics superoperator, received from context function
%
% Outputs:
%
%    fid  - a structure with four fields: fid.pos_pos, fid.pos_neg,
%           fid.neg_pos, fid.neg_neg that are used in the subsequ-
%           ent States quadrature processing
%
% Notes: the sequence starts with a pure Lz on protons at the mo-
%        ment and assumes that the relaxation superoperator is not
%        thermalised - the relaxation destination is the zero state.
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=noesyhsqc.m>

function fid=noesyhsqc(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Decouple carbon
if ismember('13C',spin_system.comp.isotopes)
    L=decouple(spin_system,L,[],{'13C'});
    R=decouple(spin_system,R,[],{'13C'});
end

% Evolution time discretization
t1.nsteps=parameters.npoints(1); t1.timestep=1/parameters.sweep(1);
t2.nsteps=parameters.npoints(2); t2.timestep=1/parameters.sweep(2);
t3.nsteps=parameters.npoints(3); t3.timestep=1/parameters.sweep(3);

% J-coupling transfer time
delta=abs(1/(4*parameters.J));

% Initial and detection states
rho=state(spin_system,'Lz','1H','cheap');
coil=state(spin_system,'L+','1H','cheap');

% Pulse operators
Np=operator(spin_system,'L+','15N'); 
Hp=operator(spin_system,'L+','1H');  
Hx=(Hp+Hp')/2; Hy=(Hp-Hp')/2i;
Nx=(Np+Np')/2; 

% Decouple nitrogen during NOESY
L_dec=decouple(spin_system,L,[],{'15N'});
R_dec=decouple(spin_system,R,[],{'15N'});

% Run NOESY forward
rho=step(spin_system,Hx,rho,pi/2);
rho_stack=evolution(spin_system,L_dec,[],rho,t1.timestep,t1.nsteps-1,'trajectory');
rho_stack_pos=coherence(spin_system,rho_stack,{{'1H',+1}});
rho_stack_neg=coherence(spin_system,rho_stack,{{'1H',-1}});
rho_stack_pos=step(spin_system,Hx,rho_stack_pos,pi/2);
rho_stack_neg=step(spin_system,Hx,rho_stack_neg,pi/2);
rho_stack_pos=homospoil(spin_system,rho_stack_pos,'destroy');
rho_stack_neg=homospoil(spin_system,rho_stack_neg,'destroy');
rho_stack_pos=evolution(spin_system,1i*R_dec,[],rho_stack_pos,parameters.tmix,1,'final');
rho_stack_neg=evolution(spin_system,1i*R_dec,[],rho_stack_neg,parameters.tmix,1,'final');
rho_stack_pos=homospoil(spin_system,rho_stack_pos,'destroy');
rho_stack_neg=homospoil(spin_system,rho_stack_neg,'destroy');
rho_stack_pos=step(spin_system,Hx,rho_stack_pos,pi/2); 
rho_stack_neg=step(spin_system,Hx,rho_stack_neg,pi/2);

% Run HSQC preparation period forward
rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,delta,1,'final');
rho_stack_neg=evolution(spin_system,L,[],rho_stack_neg,delta,1,'final');
rho_stack_pos=step(spin_system,Nx+Hx,rho_stack_pos,pi);
rho_stack_neg=step(spin_system,Nx+Hx,rho_stack_neg,pi);
rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,delta,1,'final');
rho_stack_neg=evolution(spin_system,L,[],rho_stack_neg,delta,1,'final');
rho_stack_pos=step(spin_system,Hy,rho_stack_pos,pi/2);
rho_stack_neg=step(spin_system,Hy,rho_stack_neg,pi/2);
rho_stack_pos=step(spin_system,Nx,rho_stack_pos,pi/2);
rho_stack_neg=step(spin_system,Nx,rho_stack_neg,pi/2);

% Run the rest of HSQC backward
L_dec=decouple(spin_system,L,[],{'15N'});
coil_stack=evolution(spin_system,L_dec',[],coil,-t3.timestep,t3.nsteps-1,'trajectory');
coil_stack=evolution(spin_system,L',[],coil_stack,-delta,1,'final');
coil_stack=step(spin_system,Hx+Nx,coil_stack,-pi);
coil_stack=evolution(spin_system,L',[],coil_stack,-delta,1,'final');
coil_stack=step(spin_system,Hx+Nx,coil_stack,-pi/2);
coil_stack=coherence(spin_system,coil_stack,{{'1H',0}});
coil_stack_pos=coherence(spin_system,coil_stack,{{'15N',+1}});
coil_stack_neg=coherence(spin_system,coil_stack,{{'15N',-1}});

% Stitch the halves
report(spin_system,'stitching forward and backward trajectories...');
fid.pos_pos=stitch(spin_system,L,rho_stack_pos,coil_stack_pos,{Hx},{pi},t1,t2,t3);
fid.pos_neg=stitch(spin_system,L,rho_stack_pos,coil_stack_neg,{Hx},{pi},t1,t2,t3);
fid.neg_pos=stitch(spin_system,L,rho_stack_neg,coil_stack_pos,{Hx},{pi},t1,t2,t3);
fid.neg_neg=stitch(spin_system,L,rho_stack_neg,coil_stack_neg,{Hx},{pi},t1,t2,t3);

% Permute dimensions
fid_names=fieldnames(fid);
for name_idx=1:numel(fid_names)
    fid.(fid_names{name_idx})=permute(fid.(fid_names{name_idx}),[3 2 1]);
end

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'npoints')
    error('number of points must be specificed in parameters.npoints variable.');
end
if (~isnumeric(parameters.npoints))||(~isvector(parameters.npoints))||...
   (~isreal(parameters.npoints))||(numel(parameters.npoints)~=3)||...
    any(parameters.npoints<1)||any(mod(parameters.npoints,1)~=0)
    error('parameters.npoints must be a vector of three positive integers.');
end
if ~isfield(parameters,'sweep')
    error('sweep widths must be specificed in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(~isvector(parameters.sweep))||...
   (~isreal(parameters.sweep))||(numel(parameters.sweep)~=3)
    error('parameters.sweep must be a vector of three real numbers.');
end
if ~isfield(parameters,'J')
    error('HSQC stage J-coupling must be specificed in parameters.J variable.');
end
if (~isnumeric(parameters.J))||(~isscalar(parameters.J))||...
   (~isreal(parameters.J))
    error('parameters.J must be a real number.');
end
if ~isfield(parameters,'tmix')
    error('NOESY stage mixing time must be specificed in parameters.tmix variable.');
end
if (~isnumeric(parameters.tmix))||(~isscalar(parameters.tmix))||...
   (~isreal(parameters.tmix))
    error('parameters.tmix must be a real number.');
end
end

% "The invention of the calculus of quaternions is a step towards the knowledge of
% quantities related to space which can only be compared, for its importance, with
% the invention of triple coordinates by Descartes. The ideas of this calculus, as 
% distinguished from its operations and symbols, are fitted to be of the greatest 
% use in all parts of science." - James Clerk Maxwell, 1869.
% 
% "Quaternions came from Hamilton after his really good work had been done; and,
% though beautifully ingenious, have been an unmixed evil to those who have tou-
% ched them in any way, including Clerk Maxwell." - Lord Kelvin, 1892.

