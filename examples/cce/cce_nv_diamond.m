% NV Hahn echo in diamond: the inset of Figure 3(a) in
% https://doi.org/10.1002/adts.202100254. The 493 carbon sites are the
% seed-5 realisation of the PyCCE diamond ensemble generator (40 Angstrom
% bath, [111] field); interacting pairs are retained within 6 Angstrom.
%
%                  [answer,diagnostics]=cce_nv_diamond()
%
% answer contains the CCE2 echo at each total time. diagnostics contains
% times (seconds), the CCE1 echo, and the exact isolated-nuclei product.
%
% ilya.kuprov@weizmann.ac.il

function [answer,diagnostics]=cce_nv_diamond()

% Read the occupied carbon sites, in Angstrom, with the NV axis along z
directory=fileparts(mfilename('fullpath'));
xyz=readmatrix(fullfile(directory,'nv_diamond_c13.txt'),'CommentStyle','#');
count=size(xyz,1);
gamma_c=6.7282853242747445;
gamma_e=-17608.59705;

% Match the PyCCE Zeeman convention at 500 Gauss
sys.magnet=0.05;
sys.isotopes=[{'E3'} repmat({'13C'},1,count)];
sys.output='hush';
sys.parallel={'processes',1};
shift=1e6*(gamma_c*1e7/spin('13C')-1);
inter.zeeman.matrix=[{zeros(3)} repmat({shift*eye(3)},1,count)];

% Project each electron-carbon dipolar tensor onto the NV quantisation axis
radius=vecnorm(xyz,2,2);
unit=xyz./radius;
dipolar=1.05457172*gamma_c*gamma_e*1000./(2*pi*radius.^3);
hyperfine=dipolar.*([0 0 1]-3*unit(:,3).*unit);
inter.coupling.matrix=cell(count+1);
for n=1:count
    inter.coupling.matrix{1,n+1}=[zeros(3,2) hyperfine(n,:).'];
end

% The cluster graph contains precisely the retained carbon-carbon dipoles
adjacency=false(count);
for n=1:count

    % Get neighbours within the pair cutoff, without counting pairs twice
    separation=xyz-xyz(n,:);
    distance=vecnorm(separation,2,2);
    for k=find((distance<=6)&((1:count).'>n)).'

        % Cartesian carbon-carbon coupling in Hz
        direction=separation(k,:)/distance(k);
        strength=1.05457172*gamma_c^2*1000/(2*pi*distance(k)^3);
        inter.coupling.matrix{n+1,k+1}=strength*(eye(3)-3*(direction.'*direction));
        adjacency(n,k)=true;
        adjacency(k,n)=true;

    end

end

% Native projected CCE for the 0/-1 transition
bas.formalism='zeeman-hilb';
bas.approximation='none';
parameters.observed=1;
parameters.cce.bath_spins=2:(count+1);
parameters.transition.spins=1;
parameters.transition.ms=[0 -1];
parameters.times=linspace(0,1e-4,501);
parameters.cce.order=2;
parameters.cce.approximation='projected';
parameters.cce.graph=struct('mode','explicit','adjacency',sparse(adjacency));
parameters.cce.bath_state.mode='mixed_product';
parameters.cce.basis=bas;
parameters.cce.numerics=struct('atol',1e-3,'rtol',1e-3,...
                               'max_memory_bytes',2^31);
spin_system=create(sys,inter);
[answer,details]=clexpn(spin_system,@hahn_echo,parameters,'labframe');

% Evaluate the isolated-nuclei expression independently (Eq. 29)
omega=-gamma_c*500*1000;
perp_sq=sum((2*pi*hyperfine(:,1:2)).^2,2);
conditioned=sqrt((omega-2*pi*hyperfine(:,3)).^2+perp_sq);
quarter=parameters.times/4;
isolated=prod(1-2*(perp_sq./conditioned.^2).*sin(conditioned*quarter).^2.*...
                sin(omega*quarter).^2,1).';
diagnostics=struct('times',parameters.times(:),'isolated',isolated,...
                   'bound',details.reduction.abs_error_bound,...
                   'cce1',details.reduction.by_order(:,:,2));

% Overlay the exact product with the alternating CCE1 and CCE2 markers
kfigure();
plot(1e3*diagnostics.times,isolated,'k-',...
     1e3*diagnostics.times(2:2:end),real(diagnostics.cce1(2:2:end)),'b.',...
     1e3*diagnostics.times(1:2:end),real(answer(1:2:end)),'r.');
kxlabel('total echo time / ms');
kylabel('coherence');
klegend({'isolated nuclei','CCE1','CCE2'});

end


