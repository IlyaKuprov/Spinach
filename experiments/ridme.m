% RIDME pulse sequence. Idealized hard pulses are used, the pulses only
% affect the user-specified electron. Syntax:
%
%            answer=ridme(spin_system,parameters,H,R,K)
%
% where H is the Hamiltonian commutation superoperator, R is the relaxa-
% tion superoperator and K is the chemical kinetics superoperator.
%
% Parameters:
%
%   H                         Hamiltonian (received from the context
%                             function)
%
%   R                         relaxation superoperator (received from
%                             the context function)
%
%   K                         kinetics superoperator (received from
%                             the context function)
%
%   parameters.rho0           initial state
%
%   parameters.probe_spin     number of the spin on which the
%                             sequence operates
%
%   parameters.stepsize       step size for the increment of 
%                             the relaxation period, seconds
%
%   parameters.nsteps(1)      number of steps for tau 1
% 
%   parameters.nsteps(2)      number of steps for tau 2
%
%   parameters.tmix           mixing time, seconds
%
% Outputs:
%
%   answer.pxpxpx.(real,imag)
%   answer.pypypx.(real,imag)
%   answer.mxmxpx.(real,imag)
%   answer.mymypx.(real,imag) - quadrature components of the signal
%                               corresponding to the phase cycle in-
%                               stances on third, fourth, and fifth
%                               pulse in the RIDME sequence
%
% Notes: for this experiment to work, relaxation must be present.  
%
% alice.bowen@chem.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ridme.m>

function answer=ridme(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Get operators and states
coil_imag=(state(spin_system,{'L+'},{parameters.probe_spin})+...
           state(spin_system,{'L-'},{parameters.probe_spin}))/2;
coil_real=(state(spin_system,{'L+'},{parameters.probe_spin})-...
           state(spin_system,{'L-'},{parameters.probe_spin}))/2i;
Sx=(operator(spin_system,{'L+'},{parameters.probe_spin})+...
    operator(spin_system,{'L-'},{parameters.probe_spin}))/2; 
Sy=(operator(spin_system,{'L+'},{parameters.probe_spin})-...
    operator(spin_system,{'L-'},{parameters.probe_spin}))/2i;

% First pulse
rho=step(spin_system,Sx,parameters.rho0,+pi/2);

% Evolution (tau 1)
rho_stack=evolution(spin_system,L,[],rho,(parameters.stepsize*parameters.nsteps(1)),1,'final');

% Second pulse
rho_stack=step(spin_system,Sx,rho_stack,+pi);

% Evolution (initial tau 1 and tau 2)
rho_stack=evolution(spin_system,L,[],rho_stack,parameters.stepsize,...
                   (parameters.nsteps(1)+parameters.nsteps(2)),'trajectory');

% Third pulse and a phase cycle
rho_stack_px=step(spin_system,Sx,rho_stack,+pi/2);
rho_stack_py=step(spin_system,Sy,rho_stack,+pi/2);
rho_stack_mx=step(spin_system,Sx,rho_stack,-pi/2);
rho_stack_my=step(spin_system,Sy,rho_stack,-pi/2);

% Evolution (mixing time)
rho_stack_px=evolution(spin_system,L,[],rho_stack_px,parameters.tmix,1,'final');
rho_stack_py=evolution(spin_system,L,[],rho_stack_py,parameters.tmix,1,'final');
rho_stack_mx=evolution(spin_system,L,[],rho_stack_mx,parameters.tmix,1,'final');
rho_stack_my=evolution(spin_system,L,[],rho_stack_my,parameters.tmix,1,'final');

% Fourth pulse and a phase cycle
rho_stack_pxpx=step(spin_system,Sx,rho_stack_px,+pi/2);
rho_stack_pypy=step(spin_system,Sy,rho_stack_py,+pi/2);
rho_stack_mxmx=step(spin_system,Sx,rho_stack_mx,-pi/2);
rho_stack_mymy=step(spin_system,Sy,rho_stack_my,-pi/2);

% Evolution (remaining tau 1 + tau 2)
rho_stack_pxpx(:,end:-1:1)=evolution(spin_system,L,[],rho_stack_pxpx(:,end:-1:1),parameters.stepsize,...
                                    (parameters.nsteps(1)+parameters.nsteps(2)),'refocus');
rho_stack_pypy(:,end:-1:1)=evolution(spin_system,L,[],rho_stack_pypy(:,end:-1:1),parameters.stepsize,...
                                    (parameters.nsteps(1)+parameters.nsteps(2)),'refocus');
rho_stack_mxmx(:,end:-1:1)=evolution(spin_system,L,[],rho_stack_mxmx(:,end:-1:1),parameters.stepsize,...
                                    (parameters.nsteps(1)+parameters.nsteps(2)),'refocus');
rho_stack_mymy(:,end:-1:1)=evolution(spin_system,L,[],rho_stack_mymy(:,end:-1:1),parameters.stepsize,...
                                    (parameters.nsteps(1)+parameters.nsteps(2)),'refocus');

% Fifth pulse
rho_stack_pxpxpx=step(spin_system,Sx,rho_stack_pxpx,+pi);
rho_stack_pypypx=step(spin_system,Sx,rho_stack_pypy,+pi);
rho_stack_mxmxpx=step(spin_system,Sx,rho_stack_mxmx,+pi);
rho_stack_mymypx=step(spin_system,Sx,rho_stack_mymy,+pi);

% Evolution (tau 2)
rho_stack_pxpxpx=evolution(spin_system,L,[],rho_stack_pxpxpx,(parameters.stepsize*parameters.nsteps(2)),1,'final');
rho_stack_pypypx=evolution(spin_system,L,[],rho_stack_pypypx,(parameters.stepsize*parameters.nsteps(2)),1,'final');
rho_stack_mxmxpx=evolution(spin_system,L,[],rho_stack_mxmxpx,(parameters.stepsize*parameters.nsteps(2)),1,'final');
rho_stack_mymypx=evolution(spin_system,L,[],rho_stack_mymypx,(parameters.stepsize*parameters.nsteps(2)),1,'final');

% Observation
answer.pxpxpx.real=coil_real'*rho_stack_pxpxpx/norm(coil_real,2);
answer.pxpxpx.imag=coil_imag'*rho_stack_pxpxpx/norm(coil_imag,2);
answer.pypypx.real=coil_real'*rho_stack_pypypx/norm(coil_real,2);
answer.pypypx.imag=coil_imag'*rho_stack_pypypx/norm(coil_imag,2);
answer.mxmxpx.real=coil_real'*rho_stack_mxmxpx/norm(coil_real,2);
answer.mxmxpx.imag=coil_imag'*rho_stack_mxmxpx/norm(coil_imag,2);
answer.mymypx.real=coil_real'*rho_stack_mymypx/norm(coil_real,2);
answer.mymypx.imag=coil_imag'*rho_stack_mymypx/norm(coil_imag,2);

end

% Consistency enforcement
function grumble(parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'rho0')
    error('the initial state must be provided in parameters.rho0 variable.');
end
end

% "Begin at the beginning," the King said, very gravely, 
% "and go on till you come to the end: then stop."
%
% Lewis Carroll, Alice in Wonderland

