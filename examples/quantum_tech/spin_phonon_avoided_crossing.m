% Avoided crossing between an electron spin transition and a
% quantised phonon mode in the resonant spin-phonon exchange
% model. The phonon is requested with the V# particle syntax.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function spin_phonon_avoided_crossing()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'E','V3'};

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,[]);
spin_system=basis(spin_system,bas);

% Spin operator
electron_z=operator(spin_system,'Lz',1);

% Coupling parameters
g=2*pi*4e6;
detuning=2*pi*linspace(-20e6,20e6,121);

% Exchange Hamiltonian
action=g*(operator(spin_system,{'L+','A'},{1,2})+...
          operator(spin_system,{'L-','C'},{1,2}));

% Clean up numerical asymmetry
action=(action+action')/2;

% Locate the one-excitation manifold
spin_exc=state(spin_system,{'ZL2','BL1'},{1,2});
phon_exc=state(spin_system,{'ZL1','BL2'},{1,2});
one_quant=speye(size(action,1));
spin_idx=find(diag(spin_exc)>0.5);
phon_idx=find(diag(phon_exc)>0.5);
one_quant=one_quant(:,[spin_idx phon_idx]);

% Preallocate the eigenvalue array
levels=zeros(2,numel(detuning));

% Sweep the spin-phonon detuning
for n=1:numel(detuning)

    % Build the rotating-frame Hamiltonian
    H=detuning(n)*electron_z+action;

    % Extract the dressed one-quantum doublet
    H=one_quant'*H*one_quant;
    levels(:,n)=sort(real(eig(full(H))))/(2*pi*1e6);

end

% Validate the resonant avoided-crossing gap
gap=min(levels(2,:)-levels(1,:));
if abs(gap-2*g/(2*pi*1e6))>1e-6
    error('spin-phonon avoided-crossing gap is inconsistent.');
end

% Plot the avoided crossing
kfigure(); plot(detuning/(2*pi*1e6),levels,'LineWidth',1.5);
axis tight; kgrid; kxlabel('spin-phonon detuning, MHz');
kylabel('dressed energy, MHz');
ktitle('spin-phonon avoided crossing');

end
