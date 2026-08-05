% Cross-resonance gate mechanism between two fixed-frequency
% transmons in the laboratory frame. The control transmon is
% driven at the frequency of the target transmon; through the
% static quadrature coupling, the target then rotates about an
% axis in its equatorial plane by an angle conditioned on the
% state of the control. The difference between the two condi-
% tional rotations is the ZX interaction behind the native two-
% qubit gate of fixed-frequency superconducting processors; at
% these parameters the unconditional part is the larger of the
% two, and an echo sequence would be needed to remove it. A weak
% direct tone on the target sets the gate phases. The simulation
% runs in the laboratory frame, but the recorded kets are rotated
% into the frame of the drive before the Bloch components are
% evaluated; laboratory frame transverse components would alias
% at any reasonable recording interval. Both tones run at the local
% oscillator frequency of 6.0002 GHz, which is 200 kHz above the de-
% clared bare frequency of the target transmon; that offset is the
% calibrated cross-resonance drive setting of the reference imple-
% mentation, where it is held fixed rather than derived or optimi-
% sed. Model and parameters from the cross-resonance example of the
% paraqeet package.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function cross_resonance()

% Magnet field
sys.magnet=0;

% Two transmons with three levels each
sys.isotopes={'T3','T3'};

% Transmon frequencies and anharmonicities
inter.modes.frqs={5.5e9 6.0e9};
inter.modes.anharms={-240e6 -200e6};

% Quadrature coupling between the transmons
inter.modes.exchange=cell(2,2);
inter.modes.exchange{1,2}=25e6;

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Laboratory frame drift Hamiltonian
H0=full(hamiltonian(assume(spin_system,'labframe')));

% Quadrature drive operators of the two transmons
D1=full(operator(spin_system,'C',1)+operator(spin_system,'A',1));
D2=full(operator(spin_system,'C',2)+operator(spin_system,'A',2));

% Drive tone parameters, both at the calibrated local oscillator frequency
w_lo=2*pi*6.0002e9;
amp_one=5.969026041820607e8; phi_one=0;
amp_two=2.883982055995430e7; phi_two=-0.39720756;

% Flat-top Gaussian envelope parameters
t_up=30e-9; t_down=120e-9; t_ramp=15e-9;

% Pauli observables of the target transmon qubit subspace
Sx=kron(eye(3),[0 1 0; 1 0 0; 0 0 0]);
Sy=kron(eye(3),[0 -1i 0; 1i 0 0; 0 0 0]);
Sz=kron(eye(3),[1 0 0; 0 -1 0; 0 0 0]);

% Control transmon in the ground and first excited state
unit_mat=eye(9); psi_zero=unit_mat(:,1); psi_one=unit_mat(:,4);

% Excitation number of the target transmon in the product basis
n_targ=kron(ones(3,1),(0:2)');

% Internal time stepping of the source calculation
nsub=15000; dts=1e-11; nrec=75;

% Preallocate conditional Bloch trajectories
bloch=zeros(3,nsub/nrec+1,2);
bloch(:,1,1)=[0; 0; 1]; bloch(:,1,2)=[0; 0; 1];

% Propagate on the internal midpoint grid
for k=1:nsub

    % Flat-top Gaussian envelope at the midpoint
    tmid=(k-0.5)*dts;
    env=(1+erf((tmid-t_up)/t_ramp))*(1+erf((t_down-tmid)/t_ramp))/4;

    % IQ mixer signals of the two tones
    sig_one=amp_one*env*cos(w_lo*tmid-phi_one);
    sig_two=amp_two*env*cos(w_lo*tmid-phi_two);

    % Propagate through the slice
    U=expm(-1i*(H0+sig_one*D1+sig_two*D2)*dts);
    psi_zero=U*psi_zero; psi_one=U*psi_one;

    % Record the conditional Bloch vectors of the target in the drive frame
    if mod(k,nrec)==0
        phases=exp(1i*w_lo*k*dts*n_targ);
        rot_zero=phases.*psi_zero; rot_one=phases.*psi_one;
        bloch(:,k/nrec+1,1)=real([rot_zero'*Sx*rot_zero;
                                  rot_zero'*Sy*rot_zero;
                                  rot_zero'*Sz*rot_zero]);
        bloch(:,k/nrec+1,2)=real([rot_one'*Sx*rot_one;
                                  rot_one'*Sy*rot_one;
                                  rot_one'*Sz*rot_one]);
    end

end

% Validate the norm conservation
if (abs(norm(psi_zero)-1)>1e-6)||(abs(norm(psi_one)-1)>1e-6)
    error('propagation did not conserve the norm.');
end

% Validate the signed sense of the rotation with the control in the ground state
if ~((bloch(2,end,1)>0.5)&&(bloch(3,end,1)>0.2))
    error('target rotation under the ground state control has the wrong sense.');
end

% Validate the signed sense of the rotation with the control excited
if ~((bloch(2,end,2)>0)&&(bloch(3,end,2)<-0.5))
    error('target rotation under the excited control has the wrong sense.');
end

% Plot the conditional Bloch trajectories of the target
time_axis=1e9*dts*nrec*(0:(nsub/nrec));
kfigure(); scale_figure([2.0 0.75]);
subplot(1,2,1); plot(time_axis,squeeze(bloch(:,:,1))','LineWidth',1.5);
axis tight; kgrid; kxlabel('time, ns');
kylabel('target Bloch vector, drive frame');
ktitle('control transmon in $|0\rangle$');
klegend({'$\sigma_x$','$\sigma_y$','$\sigma_z$'},'Location','Best');
subplot(1,2,2); plot(time_axis,squeeze(bloch(:,:,2))','LineWidth',1.5);
axis tight; kgrid; kxlabel('time, ns');
kylabel('target Bloch vector, drive frame');
ktitle('control transmon in $|1\rangle$');
klegend({'$\sigma_x$','$\sigma_y$','$\sigma_z$'},'Location','Best');

end

