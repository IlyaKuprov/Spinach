% Projected NV Hahn echo for the pinned PyCCE four-nucleus numerical fixture.
% Source: MICCoMpy/PyCCE d74d3cb683bf3caa07231654314e554f440833b4,
% https://doi.org/10.1002/adts.202100254. Explicit numerical fixture arrays
% are public-interface inputs, not copied implementation source. Syntax:
%
%                         [answer,diagnostics]=cce_nv_diamond()
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_nv_diamond()

% Standard Spinach system with source gyromagnetic ratios mapped to g/ppm
sys.magnet=500*1e-4; sys.isotopes={'E3','13C','13C','13C','13C'};
sys.output='hush'; sys.parallel={'local',1};
inter.zeeman.matrix=cell(1,5); inter.coupling.matrix=cell(5);
gfree=2.00231930436256;
inter.zeeman.matrix{1}=gfree*(-17608.59705*1e7/spin('E3'))*eye(3);
for n=2:5

    % Source rad/(ms Gauss) to Spinach nuclear chemical shift in ppm
    inter.zeeman.matrix{n}=1e6*(6.7282853242747445*1e7/spin('13C')-1)*eye(3);
end

% Source zero-field splitting and explicit hyperfine tensors, converted kHz to Hz
inter.coupling.matrix{1,1}=1000*diag([-960000 -960000 1920000]);
inter.coupling.matrix{1,2}=1000*[62.09198322383363 42.918773452863142 86.078663610517651;42.918773452863142 -52.625936690744069 28.63925630406008;86.078663610517651 28.63925630406008 -9.4660465330895676];
inter.coupling.matrix{1,3}=1000*[46.846946007935479 41.7131493536275 20.798152899007547;41.7131493536275 -15.576682609001635 10.413661409461005;20.798152899007547 10.413661409461005 -31.270263398933857];
inter.coupling.matrix{1,4}=1000*[7.7040985952493219 -41.664463528969883 -27.802268186621028;-41.664463528969883 7.7040985952493219 27.802268186621028;-27.802268186621028 27.802268186621028 -15.408197190498676];
inter.coupling.matrix{1,5}=1000*[14.689630582202209 -18.600579331667948 -27.874817765944968;-18.600579331667948 -13.146092067661488 13.956956441013119;-27.874817765944968 13.956956441013119 -1.5435385145407345];

% Declare the reference's observed-Sz projected approximation before propagation
for n=2:5, inter.coupling.matrix{1,n}(:,1:2)=0; end

% Explicit bath coordinates and point-dipole bath couplings in Hz
xyz=[5.3499999999999996 1.78 3.5699999999999998;7.1299999999999999 3.5699999999999998 1.78;-5.3499999999999996 5.3499999999999996 3.5699999999999998;-7.1299999999999999 3.5699999999999998 5.3499999999999996];
for j=1:4
    for k=j+1:4

        % Independent source-unit conversion; displacement is in Angstrom
        r=(xyz(k,:)-xyz(j,:))*1e-10; u=r/norm(r);
        angular=1e-7*1.054571817e-34*(6.7282853242747445*1e7)^2/norm(r)^3;
        inter.coupling.matrix{j+1,k+1}=angular*(eye(3)-3*(u.'*u))/(2*pi);
    end
end

% Local Hilbert basis for retained clusters
bas.formalism='zeeman-hilb'; bas.approximation='none';
spin_system=create(sys,inter);

% Projected transition and total echo times in seconds
parameters.observed=1; parameters.cce.bath_spins=2:5;
parameters.transition.spins=1; parameters.transition.ms=[-1 1];
parameters.times=1e-3*[0 0.012500000000000001 0.025000000000000001 0.037500000000000006 0.050000000000000003 0.0625 0.075000000000000011 0.087500000000000008 0.10000000000000001 0.1125 0.125 0.13750000000000001 0.15000000000000002 0.16250000000000001 0.17500000000000002 0.1875 0.20000000000000001 0.21250000000000002 0.22500000000000001 0.23750000000000002 0.25 0.26250000000000001 0.27500000000000002 0.28750000000000003 0.30000000000000004 0.3125 0.32500000000000001 0.33750000000000002 0.35000000000000003 0.36250000000000004 0.375 0.38750000000000001 0.40000000000000002 0.41250000000000003 0.42500000000000004 0.4375 0.45000000000000001 0.46250000000000002 0.47500000000000003 0.48750000000000004 0.5];
parameters.cce.order=2; parameters.cce.approximation='projected';
parameters.cce.graph=struct('mode','complete');
parameters.cce.bath_state.mode='mixed_product';
parameters.cce.basis=bas;
parameters.cce.numerics=struct('atol',1e-3,'rtol',1e-3,'max_memory_bytes',2^28);

% Execute the named Hahn experiment through the cluster context
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,parameters,'labframe');

% Plot the measured complex echo amplitude
kfigure(); plot(1e3*parameters.times,real(answer),'b-',...
                1e3*parameters.times,imag(answer),'r--');
xlabel('total echo time / ms'); ylabel('normalised echo amplitude E(2τ)');
legend('real','imaginary'); kgrid();

end


