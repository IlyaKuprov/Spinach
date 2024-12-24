% A MAS DNP simulation performed as described in Fred Mentink-
% Vigier's paper (Spinach rotation conventions are different):
%
%         http://dx.doi.org/10.1016/j.jmr.2015.07.001
%
% Energy level diagram as a function of the rotor phase.
%
% Calculation time: milliseconds
%
% ilya.kuprov@weizmann.ac.il

function cross_effect_mas_enlev()

% Magnet field
sys.magnet=9.394;

% Spin specification
sys.isotopes={'E','E','1H'};

% Interactions
inter.zeeman.eigs{1}=[2.0094  2.0060  2.0017];
inter.zeeman.euler{1}=[0.00 0.00 0.00];
inter.zeeman.eigs{2}=[2.0094  2.0060  2.0017];
inter.zeeman.euler{2}=pi*[107 108 124]/180;
inter.zeeman.eigs{3}=[0.00 0.00 0.00];
inter.zeeman.euler{3}=[0.00 0.00 0.00];
inter.coupling.eigs=cell(3,3);
inter.coupling.euler=cell(3,3);
inter.coupling.eigs{1,2}=[23.0e6 -11.5e6 -11.5e6];
inter.coupling.euler{1,2}=pi*[0.00 135 0.00]/180;
inter.coupling.eigs{1,3}=[1.5e6 -0.75e6 -0.75e6];
inter.coupling.euler{1,3}=[0.00 0.00 0.00];

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Stack generation parameters
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.rate=12.5e3;
parameters.rframes={};
parameters.orientation=pi*[320 141 80]/180;
parameters.masframe='magnet';
parameters.spins={'E','1H'};
parameters.offset=[0 0];
parameters.max_rank=200;

% Stack generation
H=rotor_stack(spin_system,parameters,'labframe');

% Stack diagonalization
energies=zeros(size(H{1},1),numel(H));
parfor n=1:numel(H)
    energies(:,n)=sort(real(eig(H{n})))/(2*pi*1e9);
end

% Plotting
figure(); scale_figure([0.75 2.0]);
time_axis=1e6*linspace(0,1,2*parameters.max_rank+1)'/parameters.rate;
subplot(3,1,3); plot(time_axis,energies(1:2,:)','b-'); axis tight;
kxlabel('time, $\mu$s'); kylabel('energy, GHz'); kgrid;
subplot(3,1,2); plot(time_axis,energies(3:6,:)','b-'); axis tight;
kxlabel('time, $\mu$s'); kylabel('energy, GHz'); kgrid;
subplot(3,1,1); plot(time_axis,energies(7:8,:)','b-'); axis tight;
kxlabel('time, $\mu$s'); kylabel('energy, GHz'); kgrid;

end

