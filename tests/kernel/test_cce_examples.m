% Smoke and independent-reference checks for the CCE application examples. Syntax:
%
%                    [result,reference_capture]=test_cce_examples()
%
% Outputs:
%
%     result - scalar Spinach regression-test result structure
%
%     reference_capture - scalar structure of independent physical arrays
%
% The optional capture retains independent physical arrays for regression.
%
% ilya.kuprov@weizmann.ac.il

function [result,reference_capture]=test_cce_examples()

% Make the CCE examples available for the duration of this test
saved_path=path;
path_cleanup=onCleanup(@() path(saved_path));
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'examples','cce'));

% Initialise the standard test record and preserve figure visibility
result=new_test_result('kernel/cce_examples','CCE application examples',...
                       'Runnable categories and independently evaluated small limits.');
reference_capture=struct();
visibility=get(groot,'defaultFigureVisible');
cleanup=onCleanup(@() set(groot,'defaultFigureVisible',visibility));
set(groot,'defaultFigureVisible','off');

% Pair stable reference keys with the compact calculations that produce them
cases={'cce_single_nucleus', @cce_single_nucleus;...
       'cce_nv_diamond', @cce_nv_four_spin;...
       'cce_sip_donor', @cce_sip_small;...
       'cce_sic_divacancy', @cce_sic_small;...
       'gcce_sic_clock', @gcce_clock_small;...
       'cce_p1_electron_bath', @p1_hahn_fixture;...
       'pcce_electron_layer', @pcce_patch_fixture;...
       'gcce_votpp', @votpp_density_case;...
       'cce_erbium_cawo4', @cce_erbium_cawo4_small;...
       'gcce_highspin_siv', @cce_siv_three_spin;...
       'cce_hbn_quadrupoles', @cce_hbn_quadrupoles;...
       'cce_polarized_bath', @cce_polarized_bath;...
       'gcce_two_centres', @fixture_cce_two;...
       'gcce_finite_drive', @gcce_finite_drive;...
       'cce_deer_bandwidth', @cce_deer_bandwidth;...
       'cce_noise_spectrum', @cce_noise_spectrum;...
       'mecce_dissipative_bath', @mecce_density_case;...
       'cce_methyl_molecular', @fixture_cce_methyl;...
       'cce_molecular_powder', @fixture_cce_powder};
for n=1:size(cases,1)

    % Execute each compact calculation and retain the standard harness fields only
    name=cases{n,1};
    started=tic;
    if ismember(name,{'cce_noise_spectrum','cce_single_nucleus','cce_nv_diamond',...
                          'cce_p1_electron_bath','pcce_electron_layer',...
                          'cce_erbium_cawo4','cce_hbn_quadrupoles',...
                          'cce_methyl_molecular','cce_sic_divacancy','cce_sip_donor'})
        [answer,diagnostics]=cases{n,2}();
        switch name
            case 'cce_noise_spectrum'

                % Evaluate the independent transverse correlation and retain its arrays
                omega=-spin('13C')*0.01;
                times=linspace(0,10*2*pi/abs(omega),401);
                reference=cos(omega*times(:))/4;
                x=[0 1;1 0]/2;
                reference_error=max(abs(answer-reference));
                snapshot=struct('omega',omega,'reference',reference,'x',x,'reference_error',reference_error);
            case 'cce_single_nucleus'

                % Start fresh physical inputs for this independent reference
                sys=struct();
                inter=struct();
                ref_params=struct();

                % Magnet and isotopes
                sys.magnet=0.01;
                sys.isotopes={'E3','13C'};
                sys.output='hush';
                sys.parallel={'local',1};

                % Nuclear Zeeman frequency and projected hyperfine components, rad/s
                omega=-spin('13C')*sys.magnet;
                apar=0.3*abs(omega);
                aper=0.6*abs(omega);

                % Standard Spinach tensors: electron g, nuclear ppm, and coupling Hz
                inter.zeeman.matrix={zeros(3),zeros(3)};
                inter.coupling.matrix=cell(2);
                inter.coupling.matrix{1,2}=[0 0 aper;0 0 0;0 0 apar]/(2*pi);

                % Compile metadata without constructing a full-bath basis
                create(sys,inter);

                % Projected signal; static explicit operators include their angle
                ref_params.observed_spins=1;
                ref_params.bath_spins=2;
                ref_params.expansion='standard';
                ref_params.method='cce';
                ref_params.quantity='signal';
                ref_params.order=1;
                ref_params.transition.alpha=[0;1;0];
                ref_params.transition.beta=[1;0;0];
                ref_params.assumptions='labframe';
                ref_params.orientation=[0 0 0];
                ref_params.interaction_cutoff=Inf;
                ref_params.meanfield='none';
                ref_params.graph=struct('mode','complete');
                ref_params.times=linspace(0,8*pi/abs(omega),101);

                ref_params.bath_state=struct('mode','mixed_product');
                ref_params.numerics=struct('atol',1e-9,'rtol',1e-9,'max_memory_bytes',2^28);

                % Run CCE and calculate an independent analytic echo
                ref_result=struct('value',answer,'times',ref_params.times(:),'valid',true(size(answer)));
                omega1=sqrt((omega+apar)^2+aper^2);
                reference=1-2*(aper/omega1)^2*sin(omega*ref_params.times(:)/4).^2.*...
                          sin(omega1*ref_params.times(:)/4).^2;
                ref_result.reference=reference;
                ref_result.reference_error=max(abs(ref_result.value-reference));
                fprintf('Single-nucleus analytic Hahn maximum error: %.3g\n',ref_result.reference_error);

                reference_error=ref_result.reference_error;
                snapshot=struct();
                snapshot.apar=apar;
                snapshot.aper=aper;
                snapshot.omega=omega;
                snapshot.omega1=omega1;
                snapshot.reference=reference;
                snapshot.reference_error=reference_error;
            case 'cce_nv_diamond'

                % Start fresh physical inputs for this independent reference
                sys=struct();
                inter=struct();
                ref_params=struct();

                % Standard Spinach system with source gyromagnetic ratios mapped to g/ppm
                sys.magnet=500*1e-4;
                sys.isotopes={'E3','13C','13C','13C','13C'};
                sys.output='hush';
                sys.parallel={'local',1};
                inter.zeeman.matrix=cell(1,5);
                inter.coupling.matrix=cell(5);
                gfree=2.00231930436256;
                inter.zeeman.matrix{1}=gfree*(-17608.59705*1e7/spin('E3'))*eye(3);
                for f=2:5

                    % Source rad/(ms Gauss) to Spinach nuclear chemical shift in ppm
                    inter.zeeman.matrix{f}=1e6*(6.7282853242747445*1e7/spin('13C')-1)*eye(3);

                end

                % Source zero-field splitting and explicit hyperfine tensors, converted kHz to Hz
                inter.coupling.matrix{1,1}=1000*diag([-960000 -960000 1920000]);
                inter.coupling.matrix{1,2}=1000*[62.09198322383363   42.918773452863142 86.078663610517651;...
                                              42.918773452863142 -52.625936690744069 28.63925630406008;...
                                              86.078663610517651  28.63925630406008  -9.4660465330895676];
                inter.coupling.matrix{1,3}=1000*[46.846946007935479  41.7131493536275    20.798152899007547;...
                                              41.7131493536275   -15.576682609001635  10.413661409461005;...
                                              20.798152899007547  10.413661409461005 -31.270263398933857];
                inter.coupling.matrix{1,4}=1000*[  7.7040985952493219 -41.664463528969883  -27.802268186621028;...
                                              -41.664463528969883    7.7040985952493219  27.802268186621028;...
                                              -27.802268186621028   27.802268186621028  -15.408197190498676];
                inter.coupling.matrix{1,5}=1000*[ 14.689630582202209 -18.600579331667948 -27.874817765944968;...
                                              -18.600579331667948 -13.146092067661488  13.956956441013119;...
                                              -27.874817765944968  13.956956441013119  -1.5435385145407345];

                % Declare the reference's observed-Sz projected approximation before propagation
                for f=2:5
                    inter.coupling.matrix{1,f}(:,1:2)=0;
                end

                % Explicit bath coordinates and point-dipole bath couplings in Hz
                xyz=[ 5.3499999999999996 1.78               3.5699999999999998;...
                      7.1299999999999999 3.5699999999999998 1.78;...
                     -5.3499999999999996 5.3499999999999996 3.5699999999999998;...
                     -7.1299999999999999 3.5699999999999998 5.3499999999999996];
                for j=1:4
                    for k=j+1:4

                        % Independent source-unit conversion; displacement is in Angstrom
                        r=(xyz(k,:)-xyz(j,:))*1e-10;
                        u=r/norm(r);
                        angular=1e-7*1.054571817e-34*(6.7282853242747445*1e7)^2/norm(r)^3;
                        inter.coupling.matrix{j+1,k+1}=angular*(eye(3)-3*(u.'*u))/(2*pi);

                    end
                end

                create(sys,inter);

                % The source selected two-level transition is mS=-1/+1; a selective static operator includes its angle
                ref_params.observed_spins=1;
                ref_params.bath_spins=2:5;
                ref_params.method='cce';
                ref_params.expansion='standard';
                ref_params.quantity='signal';
                ref_params.order=2;
                ref_params.transition.alpha=[0 0 1].';
                ref_params.transition.beta=[1 0 0].';
                ref_params.assumptions='labframe';
                ref_params.orientation=[0 0 0];
                ref_params.interaction_cutoff=Inf;
                ref_params.meanfield='none';
                ref_params.graph=struct('mode','complete');
                ref_params.times=1e-3*[0 0.012500000000000001 0.025000000000000001 0.037500000000000006 0.050000000000000003 0.0625 ...
                                       0.075000000000000011 0.087500000000000008 0.10000000000000001 0.1125 0.125 ...
                                       0.13750000000000001 0.15000000000000002 0.16250000000000001 0.17500000000000002 0.1875 ...
                                       0.20000000000000001 0.21250000000000002 0.22500000000000001 0.23750000000000002 0.25 ...
                                       0.26250000000000001 0.27500000000000002 0.28750000000000003 0.30000000000000004 0.3125 ...
                                       0.32500000000000001 0.33750000000000002 0.35000000000000003 0.36250000000000004 0.375 ...
                                       0.38750000000000001 0.40000000000000002 0.41250000000000003 0.42500000000000004 0.4375 ...
                                       0.45000000000000001 0.46250000000000002 0.47500000000000003 0.48750000000000004 0.5];
                a=ref_params.transition.alpha;
                b=ref_params.transition.beta;

                ref_params.bath_state=struct('mode','mixed_product');
                ref_params.numerics=struct('atol',1e-3,'rtol',1e-3,'max_memory_bytes',2^28);

                % The projected PyCCE signal uses the same fixed physical receiver
                ref_result=struct('value',answer,'times',ref_params.times(:),'valid',true(size(answer)));
                reference=[1.0000000000000073 0.80412187352127196 0.97076472650584211 0.97345787001549855 0.92208749592281747 ...
                           0.89487943408951576 0.97304683857265373 0.93078025774994988 0.86125671741742393 0.99661776302799532 ...
                           0.76355437703758766 0.99904105154528788 0.9882134199268271 0.95383133371749873 0.86766134502112169 ...
                           0.98965017367005748 0.89372993555346014 0.8732162021890324 0.98602797149751487 0.82356170821295749 ...
                           0.99419409116745483 0.93938769167121006 0.9940891681376085 0.86777114512359543 0.9860547697233486 ...
                           0.84220648442521762 0.93028570851733106 0.98826030630667416 0.82912496114994416 0.97421567083728422 ...
                           0.95985232702599865 0.98902742105821617 0.86305059097592329 0.99445386679999404 0.79888714808310557 ...
                           0.97200285714316348 0.96135900760617921 0.85148111754718481 0.96346493820569012 0.97793634578837296 ...
                           0.9615196189824301].'+...
                          1i*[3.4800759960904551e-33 -2.3669791480982818e-06 3.9053578383864199e-08 4.1251546477359897e-06 ...
                              1.3322998842642245e-06 -1.9800296667105094e-05 2.7234071561477375e-06 -4.4544694809084375e-06 ...
                              -4.2527858745940881e-06 -4.0487554818223102e-08 3.5977740510985286e-05 1.2251660190342334e-07 ...
                              6.6845336712676628e-07 8.2300776394647603e-06 -1.6303824478936545e-05 9.7554193182357096e-07 ...
                              -3.5152718524314725e-05 -4.581728372659567e-05 3.2348614343623995e-07 8.2517282974561609e-06 ...
                              -2.5910341445473197e-06 2.4206767247835291e-05 3.6414673240542315e-06 1.0680014272253585e-05 ...
                              5.5858781900956123e-06 4.8774427328959467e-05 1.3601894267783467e-05 -5.3571546260425495e-06 ...
                              -1.8314667050526693e-05 -2.0325127606064874e-05 1.7291666493708969e-05 -2.6488761044315708e-06 ...
                              -0.00015685884541566557 4.8639994138961556e-06 0.00013334140466861785 -9.8320747276684197e-06 ...
                              3.1600752332004179e-05 0.00014601325410634766 9.4604540883105856e-06 1.1986705000717478e-05 ...
                              -3.3514391192247085e-05].';
                ref_result.reference_error=max(abs(ref_result.value-reference));
                fprintf('NV pinned PyCCE CCE2 maximum complex error: %.3g.\n',ref_result.reference_error);
                fprintf('Seo ensemble T2=0.796 ms,n=2.220 is not a four-nucleus target.\n');

                reference_error=ref_result.reference_error;
                snapshot=struct();
                snapshot.a=a;
                snapshot.angular=angular;
                snapshot.b=b;
                snapshot.gfree=gfree;
                snapshot.r=r;
                snapshot.reference=reference;
                snapshot.u=u;
                snapshot.xyz=xyz;
                snapshot.reference_error=reference_error;
            case 'cce_p1_electron_bath'

                % Start fresh physical inputs for this independent reference
                sys=struct();
                inter=struct();
                ref_params=struct();

                % Magnet and isotopes
                sys.magnet=0.01;
                sys.isotopes={'E','E','14N'};
                sys.output='hush';
                sys.parallel={'local',1};

                % Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
                inter.coordinates={[0 0 0];[35 0 0];[35 0 1]};
                inter.zeeman.matrix=repmat({zeros(3)},1,3);
                inter.coupling.matrix=cell(3);

                % Explicit illustrative hyperfine and bath coupling tensors in Hz
                inter.coupling.matrix{1,2}=diag([2000 2000 -4000]);
                inter.coupling.matrix{2,3}=diag([81 81 114])*1e6;
                inter.coupling.matrix{1,3}=[0 0 0;0 0 0;400 0 200].';

                % Declare observed-Sz projection for this projected electron-bath model
                for f=2:numel(sys.isotopes)
                    if ~isempty(inter.coupling.matrix{1,f})
                        inter.coupling.matrix{1,f}(:,1:2)=0;
                    end
                end

                % Coordinates define the toy geometry; all couplings above are explicit
                inter=rmfield(inter,'coordinates');

                % Compile physical metadata without a full-bath basis
                spin_system=create(sys,inter);

                % Projected observed-spin Hahn signal and a complete small-bath graph
                ref_params.observed_spins=1;
                ref_params.bath_spins=2:3;
                ref_params.expansion='standard';
                ref_params.method='cce';
                ref_params.quantity='signal';
                ref_params.order=2;
                ref_params.transition.alpha=[1;0];
                ref_params.transition.beta=[0;1];
                ref_params.assumptions='labframe';
                ref_params.orientation=[0 0 0];
                ref_params.interaction_cutoff=Inf;
                ref_params.meanfield='none';
                ref_params.graph=struct('mode','complete');
                ref_params.times=linspace(0,200e-6,41);

                ref_params.bath_state=struct('mode','mixed_product');
                ref_params.numerics=struct('atol',1e-3,'rtol',1e-3,'max_memory_bytes',2^28);

                % Keep the electron and its nitrogen in the same physical group
                ref_params.expansion='partition';
                ref_params.groups={[2 3]};
                ref_params.partition=struct('mode','user','max_dimension',6);
                ref_params.order=1;

                % Run the cluster expansion and report its physical input boundary
                ref_result=struct('value',answer,'times',ref_params.times(:),'valid',true(size(answer)));
                fprintf('cce_p1_electron_bath: %d/%d valid values; finite-cluster demonstration.\n',...
                        nnz(ref_result.valid),numel(ref_result.valid));
                fprintf('No P1 density/T2 claim: one explicit hyperfine-bound group.\n');

                % Independently assemble projected bath Hamiltonians from input Hz tensors
                bath=ref_params.bath_spins;
                dimensions=spin_system.comp.mults(bath);
                dimension=prod(dimensions);
                operators=cell(numel(bath),3);
                for j=1:numel(bath)

                    % Angular momentum matrices in descending magnetic-quantum-number order
                    q=(dimensions(j)-1)/2;
                    m=q:-1:-q;
                    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
                    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
                    for a=1:3

                        % Embed a physical bath operator with dense Kronecker products
                        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
                                            kron(local{a},eye(prod(dimensions(j+1:end)))));

                    end

                end
                bath_h=zeros(dimension);
                difference=zeros(dimension);
                for j=1:numel(bath)

                    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
                    id=bath(j);
                    if sys.isotopes{id}(1)~='E'
                        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
                    end
                    for a=1:3

                        % Projected hyperfine difference and on-site quadratic terms
                        if ~isempty(inter.coupling.matrix{1,id})
                            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
                        end
                        for b=1:3
                            if ~isempty(inter.coupling.matrix{id,id})
                                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
                            end
                        end

                    end
                    for k=j+1:numel(bath)
                        for a=1:3
                            for b=1:3

                                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                                tensor=inter.coupling.matrix{id,bath(k)};
                                if ~isempty(tensor)
                                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                                end

                            end
                        end
                    end

                end
                ha=bath_h+difference/2;
                hb=bath_h-difference/2;
                reference=complex(zeros(numel(ref_result.times),1));
                for f=1:numel(ref_result.times)

                    % Independent full-bath Hahn propagators for a maximally mixed state
                    ua=expm(-1i*ha*ref_result.times(f)/2);
                    ub=expm(-1i*hb*ref_result.times(f)/2);
                    reference(f)=conj(trace(ub*ua*(ua*ub)'))/dimension;

                end
                ref_result.reference_error=max(abs(ref_result.value-reference));
                fprintf('Independent dense Hahn maximum complex error: %.3g.\n',ref_result.reference_error);

                reference_error=ref_result.reference_error;
                snapshot=struct();
                snapshot.a=a;
                snapshot.b=b;
                snapshot.bath=bath;
                snapshot.bath_h=bath_h;
                snapshot.difference=difference;
                snapshot.dimension=dimension;
                snapshot.dimensions=dimensions;
                snapshot.ha=ha;
                snapshot.hb=hb;
                snapshot.id=id;
                snapshot.local=local;
                snapshot.m=m;
                snapshot.operators=operators;
                snapshot.plus=plus;
                snapshot.reference=reference;
                snapshot.tensor=tensor;
                snapshot.ua=ua;
                snapshot.ub=ub;
                snapshot.reference_error=reference_error;
            case 'pcce_electron_layer'

                % Start fresh physical inputs for this independent reference
                sys=struct();
                inter=struct();
                ref_params=struct();

                % Magnet and isotopes
                sys.magnet=0.01;
                sys.isotopes={'E','E','E','E'};
                sys.output='hush';
                sys.parallel={'local',1};

                % Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
                inter.coordinates={[0 0 0];[20 0 0];[21 0 0];[0 30 0]};
                inter.zeeman.matrix=repmat({zeros(3)},1,4);
                inter.coupling.matrix=cell(4);

                % Explicit illustrative hyperfine and bath coupling tensors in Hz
                inter.coupling.matrix{1,2}=[0 0 0;0 0 0;15000 0 8000].';
                inter.coupling.matrix{1,3}=[0 0 0;0 0 0;-12000 0 6000].';
                inter.coupling.matrix{1,4}=[0 0 0;0 0 0;4000 0 3000].';
                inter.coupling.matrix{2,3}=diag([2e5 2e5 -4e5]);
                inter.coupling.matrix{3,4}=diag([300 300 -600]);

                % Declare observed-Sz projection for this projected electron-bath model
                for f=2:numel(sys.isotopes)
                    if ~isempty(inter.coupling.matrix{1,f})
                        inter.coupling.matrix{1,f}(:,1:2)=0;
                    end
                end

                % Coordinates define the toy geometry; all couplings above are explicit
                inter=rmfield(inter,'coordinates');

                % Compile physical metadata without a full-bath basis
                spin_system=create(sys,inter);

                % Projected observed-spin Hahn signal and a complete small-bath graph
                ref_params.observed_spins=1;
                ref_params.bath_spins=2:4;
                ref_params.expansion='standard';
                ref_params.method='cce';
                ref_params.quantity='signal';
                ref_params.order=3;
                ref_params.transition.alpha=[1;0];
                ref_params.transition.beta=[0;1];
                ref_params.assumptions='labframe';
                ref_params.orientation=[0 0 0];
                ref_params.interaction_cutoff=Inf;
                ref_params.meanfield='none';
                ref_params.graph=struct('mode','complete');
                ref_params.times=linspace(0,200e-6,41);

                ref_params.bath_state=struct('mode','mixed_product');
                ref_params.numerics=struct('atol',1e-5,'rtol',1e-5,'max_memory_bytes',2^28);

                % Freeze the strongly coupled pair as one partition
                ref_params.expansion='partition';
                ref_params.groups={[2 3],4};
                ref_params.partition=struct('mode','user','max_dimension',4);
                ref_params.order=2;

                % Run the cluster expansion and report its physical input boundary
                ref_result=struct('value',answer,'times',ref_params.times(:),'valid',true(size(answer)));
                fprintf('pcce_electron_layer: %d/%d valid values; finite-cluster demonstration.\n',...
                        nnz(ref_result.valid),numel(ref_result.valid));

                % Compare with a separate unpartitioned full-order expansion
                full_parameters=rmfield(ref_params,{'partition','groups'});
                full_parameters.expansion='standard';
                full_parameters.order=3;
                native_parameters=struct('observed',1,'times',ref_params.times,...
                                         'transition',struct('spins',1,'ms',[1/2 -1/2]),...
                                         'cce',struct('bath_spins',2:4,'order',3,'approximation','projected',...
                                                      'graph',ref_params.graph,'numerics',ref_params.numerics));
                [full_result.value,full_diagnostics]=clexpn(spin_system,@hahn_echo,native_parameters,'labframe');
                ref_result.reference_error=max(abs(ref_result.value-full_result.value));
                fprintf('Partition/full-order maximum difference: %.3g (same Hamiltonian).\n',...
                        ref_result.reference_error);

                % Independently assemble projected bath Hamiltonians from input Hz tensors
                bath=ref_params.bath_spins;
                dimensions=spin_system.comp.mults(bath);
                dimension=prod(dimensions);
                operators=cell(numel(bath),3);
                for j=1:numel(bath)

                    % Angular momentum matrices in descending magnetic-quantum-number order
                    q=(dimensions(j)-1)/2;
                    m=q:-1:-q;
                    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
                    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
                    for a=1:3

                        % Embed a physical bath operator with dense Kronecker products
                        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
                                            kron(local{a},eye(prod(dimensions(j+1:end)))));

                    end

                end
                bath_h=zeros(dimension);
                difference=zeros(dimension);
                for j=1:numel(bath)

                    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
                    id=bath(j);
                    if sys.isotopes{id}(1)~='E'
                        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
                    end
                    for a=1:3

                        % Projected hyperfine difference and on-site quadratic terms
                        if ~isempty(inter.coupling.matrix{1,id})
                            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
                        end
                        for b=1:3
                            if ~isempty(inter.coupling.matrix{id,id})
                                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
                            end
                        end

                    end
                    for k=j+1:numel(bath)
                        for a=1:3
                            for b=1:3

                                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                                tensor=inter.coupling.matrix{id,bath(k)};
                                if ~isempty(tensor)
                                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                                end

                            end
                        end
                    end

                end
                ha=bath_h+difference/2;
                hb=bath_h-difference/2;
                reference=complex(zeros(numel(ref_result.times),1));
                for f=1:numel(ref_result.times)

                    % Independent full-bath Hahn propagators for a maximally mixed state
                    ua=expm(-1i*ha*ref_result.times(f)/2);
                    ub=expm(-1i*hb*ref_result.times(f)/2);
                    reference(f)=conj(trace(ub*ua*(ua*ub)'))/dimension;

                end
                ref_result.reference_error=max(abs(ref_result.value-reference));
                fprintf('Independent dense Hahn maximum complex error: %.3g.\n',ref_result.reference_error);

                reference_error=ref_result.reference_error;
                snapshot=struct();
                snapshot.a=a;
                snapshot.b=b;
                snapshot.bath=bath;
                snapshot.bath_h=bath_h;
                snapshot.difference=difference;
                snapshot.dimension=dimension;
                snapshot.dimensions=dimensions;
                snapshot.full_parameters=full_parameters;
                snapshot.ha=ha;
                snapshot.hb=hb;
                snapshot.id=id;
                snapshot.local=local;
                snapshot.m=m;
                snapshot.operators=operators;
                snapshot.plus=plus;
                snapshot.reference=reference;
                snapshot.tensor=tensor;
                snapshot.ua=ua;
                snapshot.ub=ub;
                snapshot.reference_error=reference_error;
                snapshot.full_order_error=max(abs(answer-full_result.value));
                snapshot.full_answer=full_result.value;
                snapshot.full_diagnostics=full_diagnostics;
                result=test_true(result,'Partition full-order identity',...
                                 snapshot.full_order_error<1e-9,...
                                 'The retained pair agrees with the unpartitioned full bath.');
            case 'cce_erbium_cawo4'

                % Start fresh physical inputs for this independent reference
                sys=struct();
                inter=struct();
                ref_params=struct();

                % Magnet and isotopes
                sys.magnet=0.01;
                sys.isotopes={'E','183W','183W'};
                sys.output='hush';
                sys.parallel={'local',1};

                % Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
                inter.coordinates={[0 0 0];[8 0 2];[0 10 -1]};
                inter.zeeman.matrix=repmat({zeros(3)},1,3);
                inter.coupling.matrix=cell(3);

                % Explicit illustrative hyperfine and bath coupling tensors in Hz
                inter.coupling.matrix{1,2}=[0 0 0;0 0 0;3000 1200 2000].';
                inter.coupling.matrix{1,3}=[0 0 0;0 0 0;-900 400 1000].';
                inter.coupling.matrix{2,3}=diag([80 80 -160]);

                % Coordinates define the toy geometry; all couplings above are explicit
                inter=rmfield(inter,'coordinates');

                spin_system=create(sys,inter);
                ref_params.bath=2:3;
                ref_params.times=linspace(0,200e-6,41);

                % Independently assemble projected bath Hamiltonians from input Hz tensors
                bath=ref_params.bath;
                dimensions=spin_system.comp.mults(bath);
                dimension=prod(dimensions);
                operators=cell(numel(bath),3);
                for j=1:numel(bath)

                    % Angular momentum matrices in descending magnetic-quantum-number order
                    q=(dimensions(j)-1)/2;
                    m=q:-1:-q;
                    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
                    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
                    for a=1:3

                        % Embed a physical bath operator with dense Kronecker products
                        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
                                            kron(local{a},eye(prod(dimensions(j+1:end)))));

                    end

                end
                bath_h=zeros(dimension);
                difference=zeros(dimension);
                for j=1:numel(bath)

                    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
                    id=bath(j);
                    if sys.isotopes{id}(1)~='E'
                        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
                    end
                    for a=1:3

                        % Projected hyperfine difference and on-site quadratic terms
                        if ~isempty(inter.coupling.matrix{1,id})
                            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
                        end
                        for b=1:3
                            if ~isempty(inter.coupling.matrix{id,id})
                                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
                            end
                        end

                    end
                    for k=j+1:numel(bath)
                        for a=1:3
                            for b=1:3

                                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                                tensor=inter.coupling.matrix{id,bath(k)};
                                if ~isempty(tensor)
                                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                                end

                            end
                        end
                    end

                end
                ha=bath_h+difference/2;
                hb=bath_h-difference/2;
                reference=complex(zeros(numel(ref_params.times),1));
                for f=1:numel(ref_params.times)

                    % Independent full-bath Hahn propagators for a maximally mixed state
                    ua=expm(-1i*ha*ref_params.times(f)/2);
                    ub=expm(-1i*hb*ref_params.times(f)/2);
                    reference(f)=conj(trace(ub*ua*(ua*ub)'))/dimension;

                end
                reference_error=max(abs(answer-reference));

                snapshot=struct();
                snapshot.a=a;
                snapshot.b=b;
                snapshot.bath=bath;
                snapshot.dimensions=dimensions;
                snapshot.dimension=dimension;
                snapshot.operators=operators;
                snapshot.q=q;
                snapshot.m=m;
                snapshot.plus=plus;
                snapshot.local=local;
                snapshot.bath_h=bath_h;
                snapshot.difference=difference;
                snapshot.id=id;
                snapshot.tensor=tensor;
                snapshot.ha=ha;
                snapshot.hb=hb;
                snapshot.reference=reference;
                snapshot.ua=ua;
                snapshot.ub=ub;
                snapshot.reference_error=reference_error;
            case 'cce_hbn_quadrupoles'

                % Start fresh physical inputs for this independent reference
                sys=struct();
                inter=struct();
                ref_params=struct();

                % Magnet and isotopes
                sys.magnet=0.01;
                sys.isotopes={'E','11B','14N'};
                sys.output='hush';
                sys.parallel={'local',1};

                % Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
                inter.coordinates={[0 0 0];[1.45 0 0];[-0.725 1.2557 0]};
                inter.zeeman.matrix=repmat({zeros(3)},1,3);
                inter.coupling.matrix=cell(3);

                % Explicit illustrative hyperfine and bath coupling tensors in Hz
                inter.coupling.matrix{1,2}=[0 0 0;0 0 0;40000 0 20000].';
                inter.coupling.matrix{1,3}=[0 0 0;0 0 0;20000 0 -10000].';
                inter.coupling.matrix{2,3}=diag([500 500 -1000]);
                inter.coupling.matrix{2,2}=diag([-1 -1 2])*2e5;
                inter.coupling.matrix{3,3}=diag([-1 -1 2])*5e4;

                % Coordinates define the toy geometry; all couplings above are explicit
                inter=rmfield(inter,'coordinates');

                spin_system=create(sys,inter);
                ref_params.bath=2:3;
                ref_params.times=linspace(0,200e-6,41);

                % Independently assemble projected bath Hamiltonians from input Hz tensors
                bath=ref_params.bath;
                dimensions=spin_system.comp.mults(bath);
                dimension=prod(dimensions);
                operators=cell(numel(bath),3);
                for j=1:numel(bath)

                    % Angular momentum matrices in descending magnetic-quantum-number order
                    q=(dimensions(j)-1)/2;
                    m=q:-1:-q;
                    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
                    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
                    for a=1:3

                        % Embed a physical bath operator with dense Kronecker products
                        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
                                            kron(local{a},eye(prod(dimensions(j+1:end)))));

                    end

                end
                bath_h=zeros(dimension);
                difference=zeros(dimension);
                for j=1:numel(bath)

                    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
                    id=bath(j);
                    if sys.isotopes{id}(1)~='E'
                        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
                    end
                    for a=1:3

                        % Projected hyperfine difference and on-site quadratic terms
                        if ~isempty(inter.coupling.matrix{1,id})
                            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
                        end
                        for b=1:3
                            if ~isempty(inter.coupling.matrix{id,id})
                                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
                            end
                        end

                    end
                    for k=j+1:numel(bath)
                        for a=1:3
                            for b=1:3

                                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                                tensor=inter.coupling.matrix{id,bath(k)};
                                if ~isempty(tensor)
                                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                                end

                            end
                        end
                    end

                end
                ha=bath_h+difference/2;
                hb=bath_h-difference/2;
                reference=complex(zeros(numel(ref_params.times),1));
                for f=1:numel(ref_params.times)

                    % Independent full-bath Hahn propagators for a maximally mixed state
                    ua=expm(-1i*ha*ref_params.times(f)/2);
                    ub=expm(-1i*hb*ref_params.times(f)/2);
                    reference(f)=conj(trace(ub*ua*(ua*ub)'))/dimension;

                end
                reference_error=max(abs(answer-reference));

                snapshot=struct();
                snapshot.a=a;
                snapshot.b=b;
                snapshot.bath=bath;
                snapshot.dimensions=dimensions;
                snapshot.dimension=dimension;
                snapshot.operators=operators;
                snapshot.q=q;
                snapshot.m=m;
                snapshot.plus=plus;
                snapshot.local=local;
                snapshot.bath_h=bath_h;
                snapshot.difference=difference;
                snapshot.id=id;
                snapshot.tensor=tensor;
                snapshot.ha=ha;
                snapshot.hb=hb;
                snapshot.reference=reference;
                snapshot.ua=ua;
                snapshot.ub=ub;
                snapshot.reference_error=reference_error;
            case 'cce_methyl_molecular'

                % Start fresh physical inputs for this independent reference
                sys=struct();
                inter=struct();
                ref_params=struct();

                % Magnet and isotopes
                sys.magnet=0.01;
                sys.isotopes={'E','1H','1H','1H'};
                sys.output='hush';
                sys.parallel={'local',1};

                % Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
                inter.coordinates={[0 0 0];[5 1 0];[5 -0.5 0.866];[5 -0.5 -0.866]};
                inter.zeeman.matrix=repmat({zeros(3)},1,4);
                inter.coupling.matrix=cell(4);

                % Explicit illustrative hyperfine and bath coupling tensors in Hz
                inter.coupling.matrix{1,2}=[0 0 0;0 0 0;7000 0 5000].';
                inter.coupling.matrix{1,3}=[0 0 0;0 0 0;5000 0 4000].';
                inter.coupling.matrix{1,4}=[0 0 0;0 0 0;6000 0 3000].';
                inter.coupling.matrix{2,3}=2000*eye(3);
                inter.coupling.matrix{2,4}=2000*eye(3);
                inter.coupling.matrix{3,4}=2000*eye(3);

                % Coordinates define the toy geometry; all couplings above are explicit
                inter=rmfield(inter,'coordinates');

                spin_system=create(sys,inter);
                ref_params.bath=2:4;
                ref_params.times=linspace(0,200e-6,41);

                % Independently assemble projected bath Hamiltonians from input Hz tensors
                bath=ref_params.bath;
                dimensions=spin_system.comp.mults(bath);
                dimension=prod(dimensions);
                operators=cell(numel(bath),3);
                for j=1:numel(bath)

                    % Angular momentum matrices in descending magnetic-quantum-number order
                    q=(dimensions(j)-1)/2;
                    m=q:-1:-q;
                    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
                    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
                    for a=1:3

                        % Embed a physical bath operator with dense Kronecker products
                        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
                                            kron(local{a},eye(prod(dimensions(j+1:end)))));

                    end

                end
                bath_h=zeros(dimension);
                difference=zeros(dimension);
                for j=1:numel(bath)

                    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
                    id=bath(j);
                    if sys.isotopes{id}(1)~='E'
                        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
                    end
                    for a=1:3

                        % Projected hyperfine difference and on-site quadratic terms
                        if ~isempty(inter.coupling.matrix{1,id})
                            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
                        end
                        for b=1:3
                            if ~isempty(inter.coupling.matrix{id,id})
                                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
                            end
                        end

                    end
                    for k=j+1:numel(bath)
                        for a=1:3
                            for b=1:3

                                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                                tensor=inter.coupling.matrix{id,bath(k)};
                                if ~isempty(tensor)
                                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                                end

                            end
                        end
                    end

                end
                ha=bath_h+difference/2;
                hb=bath_h-difference/2;
                reference=complex(zeros(numel(ref_params.times),1));
                for f=1:numel(ref_params.times)

                    % Independent full-bath Hahn propagators for a maximally mixed state
                    ua=expm(-1i*ha*ref_params.times(f)/2);
                    ub=expm(-1i*hb*ref_params.times(f)/2);
                    reference(f)=conj(trace(ub*ua*(ua*ub)'))/dimension;

                end
                reference_error=max(abs(answer-reference));

                snapshot=struct();
                snapshot.a=a;
                snapshot.b=b;
                snapshot.bath=bath;
                snapshot.dimensions=dimensions;
                snapshot.dimension=dimension;
                snapshot.operators=operators;
                snapshot.q=q;
                snapshot.m=m;
                snapshot.plus=plus;
                snapshot.local=local;
                snapshot.bath_h=bath_h;
                snapshot.difference=difference;
                snapshot.id=id;
                snapshot.tensor=tensor;
                snapshot.ha=ha;
                snapshot.hb=hb;
                snapshot.reference=reference;
                snapshot.ua=ua;
                snapshot.ub=ub;
                snapshot.reference_error=reference_error;
            case 'cce_sic_divacancy'

                % Start fresh physical inputs for this independent reference
                sys=struct();
                inter=struct();
                ref_params=struct();

                % Magnet and isotopes
                sys.magnet=0.01;
                sys.isotopes={'E','13C','29Si'};
                sys.output='hush';
                sys.parallel={'local',1};

                % Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
                inter.coordinates={[0 0 0];[4 0 1];[0 5 2]};
                inter.zeeman.matrix=repmat({zeros(3)},1,3);
                inter.coupling.matrix=cell(3);

                % Explicit illustrative hyperfine and bath coupling tensors in Hz
                inter.coupling.matrix{1,2}=[0 0 0;0 0 0;15000 0 8000].';
                inter.coupling.matrix{1,3}=[0 0 0;0 0 0;9000 0 -4000].';
                inter.coupling.matrix{2,3}=diag([250 250 -500]);

                % Coordinates define the toy geometry; all couplings above are explicit
                inter=rmfield(inter,'coordinates');

                spin_system=create(sys,inter);
                ref_params.bath=2:3;
                ref_params.times=linspace(0,200e-6,41);

                % Independently assemble projected bath Hamiltonians from input Hz tensors
                bath=ref_params.bath;
                dimensions=spin_system.comp.mults(bath);
                dimension=prod(dimensions);
                operators=cell(numel(bath),3);
                for j=1:numel(bath)

                    % Angular momentum matrices in descending magnetic-quantum-number order
                    q=(dimensions(j)-1)/2;
                    m=q:-1:-q;
                    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
                    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
                    for a=1:3

                        % Embed a physical bath operator with dense Kronecker products
                        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
                                            kron(local{a},eye(prod(dimensions(j+1:end)))));

                    end

                end
                bath_h=zeros(dimension);
                difference=zeros(dimension);
                for j=1:numel(bath)

                    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
                    id=bath(j);
                    if sys.isotopes{id}(1)~='E'
                        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
                    end
                    for a=1:3

                        % Projected hyperfine difference and on-site quadratic terms
                        if ~isempty(inter.coupling.matrix{1,id})
                            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
                        end
                        for b=1:3
                            if ~isempty(inter.coupling.matrix{id,id})
                                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
                            end
                        end

                    end
                    for k=j+1:numel(bath)
                        for a=1:3
                            for b=1:3

                                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                                tensor=inter.coupling.matrix{id,bath(k)};
                                if ~isempty(tensor)
                                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                                end

                            end
                        end
                    end

                end
                ha=bath_h+difference/2;
                hb=bath_h-difference/2;
                reference=complex(zeros(numel(ref_params.times),1));
                for f=1:numel(ref_params.times)

                    % Independent full-bath Hahn propagators for a maximally mixed state
                    ua=expm(-1i*ha*ref_params.times(f)/2);
                    ub=expm(-1i*hb*ref_params.times(f)/2);
                    reference(f)=conj(trace(ub*ua*(ua*ub)'))/dimension;

                end
                reference_error=max(abs(answer-reference));

                snapshot=struct();
                snapshot.a=a;
                snapshot.b=b;
                snapshot.bath=bath;
                snapshot.dimensions=dimensions;
                snapshot.dimension=dimension;
                snapshot.operators=operators;
                snapshot.q=q;
                snapshot.m=m;
                snapshot.plus=plus;
                snapshot.local=local;
                snapshot.bath_h=bath_h;
                snapshot.difference=difference;
                snapshot.id=id;
                snapshot.tensor=tensor;
                snapshot.ha=ha;
                snapshot.hb=hb;
                snapshot.reference=reference;
                snapshot.ua=ua;
                snapshot.ub=ub;
                snapshot.reference_error=reference_error;
            case 'cce_sip_donor'

                % Start fresh physical inputs for this independent reference
                sys=struct();
                inter=struct();
                ref_params=struct();

                % Magnet and isotopes
                sys.magnet=0.01;
                sys.isotopes={'E','29Si','29Si'};
                sys.output='hush';
                sys.parallel={'local',1};

                % Coordinates in Angstrom, zero electron g, nuclear shifts in ppm
                inter.coordinates={[0 0 0];[8 4 0];[12 2 4]};
                inter.zeeman.matrix=repmat({zeros(3)},1,3);
                inter.coupling.matrix=cell(3);

                % Explicit illustrative hyperfine and bath coupling tensors in Hz
                a=25;
                k0=0.85*2*pi/5.43;
                for f=2:3

                    % Hydrogenic envelope, and six-valley interference, illustrative scale
                    r=inter.coordinates{f};
                    contact=2e5*exp(-2*norm(r)/a)*(sum(cos(k0*r))/3)^2;

                    % Retain observed Sz only for the explicitly projected donor model
                    inter.coupling.matrix{1,f}=diag([0 0 contact]);

                end
                inter.coupling.matrix{2,3}=diag([200 200 -400]);

                % Coordinates define the toy geometry; all couplings above are explicit
                inter=rmfield(inter,'coordinates');

                spin_system=create(sys,inter);
                ref_params.bath=2:3;
                ref_params.times=linspace(0,200e-6,41);

                % Independently assemble projected bath Hamiltonians from input Hz tensors
                bath=ref_params.bath;
                dimensions=spin_system.comp.mults(bath);
                dimension=prod(dimensions);
                operators=cell(numel(bath),3);
                for j=1:numel(bath)

                    % Angular momentum matrices in descending magnetic-quantum-number order
                    q=(dimensions(j)-1)/2;
                    m=q:-1:-q;
                    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
                    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
                    for a=1:3

                        % Embed a physical bath operator with dense Kronecker products
                        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
                                            kron(local{a},eye(prod(dimensions(j+1:end)))));

                    end

                end
                bath_h=zeros(dimension);
                difference=zeros(dimension);
                for j=1:numel(bath)

                    % Nuclear Zeeman terms; the illustrative electron g tensors are zero
                    id=bath(j);
                    if sys.isotopes{id}(1)~='E'
                        bath_h=bath_h-spin(sys.isotopes{id})*sys.magnet*operators{j,3};
                    end
                    for a=1:3

                        % Projected hyperfine difference and on-site quadratic terms
                        if ~isempty(inter.coupling.matrix{1,id})
                            difference=difference+2*pi*inter.coupling.matrix{1,id}(a,3)*operators{j,a};
                        end
                        for b=1:3
                            if ~isempty(inter.coupling.matrix{id,id})
                                bath_h=bath_h+2*pi*inter.coupling.matrix{id,id}(a,b)*operators{j,a}*operators{j,b};
                            end
                        end

                    end
                    for k=j+1:numel(bath)
                        for a=1:3
                            for b=1:3

                                % Explicit pair tensors, with no call to the CCE Hamiltonian builder
                                tensor=inter.coupling.matrix{id,bath(k)};
                                if ~isempty(tensor)
                                    bath_h=bath_h+2*pi*tensor(a,b)*operators{j,a}*operators{k,b};
                                end

                            end
                        end
                    end

                end
                ha=bath_h+difference/2;
                hb=bath_h-difference/2;
                reference=complex(zeros(numel(ref_params.times),1));
                for f=1:numel(ref_params.times)

                    % Independent full-bath Hahn propagators for a maximally mixed state
                    ua=expm(-1i*ha*ref_params.times(f)/2);
                    ub=expm(-1i*hb*ref_params.times(f)/2);
                    reference(f)=conj(trace(ub*ua*(ua*ub)'))/dimension;

                end
                reference_error=max(abs(answer-reference));

                snapshot=struct();
                snapshot.a=a;
                snapshot.b=b;
                snapshot.bath=bath;
                snapshot.dimensions=dimensions;
                snapshot.dimension=dimension;
                snapshot.operators=operators;
                snapshot.q=q;
                snapshot.m=m;
                snapshot.plus=plus;
                snapshot.local=local;
                snapshot.bath_h=bath_h;
                snapshot.difference=difference;
                snapshot.id=id;
                snapshot.tensor=tensor;
                snapshot.ha=ha;
                snapshot.hb=hb;
                snapshot.reference=reference;
                snapshot.ua=ua;
                snapshot.ub=ub;
                snapshot.reference_error=reference_error;
                snapshot.r=r;
                snapshot.contact=contact;
                snapshot.k0=k0;
        end
        actual=struct('value',answer,'reference_error',reference_error,...
                      'valid',diagnostics.reduction.valid,'density',[]);
        reference_capture.(name)=snapshot;
        reference_capture.(name).diagnostics=diagnostics;
        result=test_true(result,[name ' projected error bounds'],...
                         all(diagnostics.reduction.abs_error_bound>0),...
                         'Projected spectral error bounds must reach the reducer.');
    elseif ismember(name,{'gcce_sic_clock','gcce_highspin_siv','gcce_votpp',...
                              'gcce_two_centres','gcce_finite_drive','mecce_dissipative_bath'})
        [answer,diagnostics]=cases{n,2}();
        actual=diagnostics.reduction;
        actual.times=diagnostics.parameters.times(:);

        % Check the observed density without repairing its populations
        dimension=sqrt(size(answer,2));
        actual.density=reshape(answer.',dimension,dimension,[]);
        actual.diagnostics.trace_error=zeros(size(answer,1),1);
        actual.diagnostics.minimum_eigenvalue=zeros(size(answer,1),1);
        for sample=1:size(answer,1)
            rho=actual.density(:,:,sample);
            actual.diagnostics.trace_error(sample)=abs(trace(rho)-1);
            actual.diagnostics.minimum_eigenvalue(sample)=min(eig((rho+rho')/2));
        end
        switch name
            case 'gcce_sic_clock'

                % Start fresh physical inputs for this independent reference
                sys=struct();
                inter=struct();

                % Magnet and isotopes; zero electron g isolates the model energy scales
                sys.magnet=0;
                sys.isotopes={'E3','13C','29Si'};
                sys.output='hush';
                sys.parallel={'local',1};
                inter.zeeman.matrix=repmat({zeros(3)},1,3);
                inter.coupling.matrix=cell(3);

                % Standard Spinach coupling and zero-field-splitting tensors in Hz
                inter.coupling.matrix{1,1}=diag([-1.0e4+2e3 -1.0e4-2e3 2e4]);
                inter.coupling.matrix{1,2}=diag([2200 1800 3000]);
                inter.coupling.matrix{1,3}=diag([900 1200 1700]);
                inter.coupling.matrix{2,3}=diag([120 120 -240]);

                % Compile physical metadata and the independently specified preparation
                spin_system=create(sys,inter);
                v=[1;0;1]/sqrt(2);
                observed=v*v';
                times=linspace(0,200e-6,31);

                % Independently construct the full dense Hamiltonian from the input Hz tensors
                dimensions=spin_system.comp.mults;
                dimension=prod(dimensions);
                operators=cell(numel(dimensions),3);
                h=zeros(dimension);
                for j=1:numel(dimensions)

                    % Local angular momentum matrices and their tensor-product embeddings
                    q=(dimensions(j)-1)/2;
                    m=q:-1:-q;
                    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
                    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
                    for a=1:3
                        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
                                            kron(local{a},eye(prod(dimensions(j+1:end)))));
                    end

                end
                for j=1:numel(dimensions)
                    for k=j:numel(dimensions)

                        % Spinach ordered pair matrices transpose into Cartesian operators
                        tensor=inter.coupling.matrix{j,k};
                        if isempty(tensor)
                            continue;
                        end
                        for a=1:3
                            for b=1:3
                                h=h+2*pi*tensor(b,a)*operators{j,a}*operators{k,b};
                            end
                        end

                    end
                end

                % Propagate an explicit product density and trace over the whole finite bath
                d_obs=size(observed,1);
                db=dimension/d_obs;
                rho0=kron(observed,eye(db)/db);
                reference=zeros(size(answer));
                for f=1:numel(times)

                    % This independent calculation uses expm, not the CCE eigensolver
                    u=expm(-1i*h*times(f));
                    rho=u*rho0*u';
                    reduced=zeros(d_obs);
                    for k=1:db
                        ids=k:db:dimension;
                        reduced=reduced+rho(ids,ids);
                    end
                    reference(f,:)=reduced(:).';

                end
                reference_error=max(abs(answer-reference),[],'all');

                % Preserve independent physical arrays for the stage-0 comparison
                snapshot=struct('a',a,'b',b,'db',db,'d_obs',d_obs,...
                                'dimension',dimension,'dimensions',dimensions,'h',h,...
                                'local',{local},'m',m,'operators',{operators},...
                                'plus',plus,'reduced',reduced,'reference',reference,...
                                'rho',rho,'rho0',rho0,'u',u,'v',v,'observed',observed,...
                                'reference_error',reference_error,'tensor',tensor,'ids',ids);

                % Retain this case reference in the public capture structure
                actual.reference_error=reference_error;
                reference_capture.gcce_sic_clock=snapshot;
            case 'gcce_highspin_siv'

                % Start fresh physical inputs for this independent reference
                sys=struct();
                inter=struct();

                % Magnet and isotopes; zero electron g isolates the model energy scales
                sys.magnet=0;
                sys.isotopes={'E4','13C','29Si'};
                sys.output='hush';
                sys.parallel={'local',1};
                inter.zeeman.matrix=repmat({zeros(3)},1,3);
                inter.coupling.matrix=cell(3);

                % Standard Spinach coupling and zero-field-splitting tensors in Hz
                inter.coupling.matrix{1,1}=diag([-1 -1 2])*1e4;
                inter.coupling.matrix{1,2}=diag([1700 1300 2400]);
                inter.coupling.matrix{1,3}=diag([800 900 1500]);
                inter.coupling.matrix{2,3}=diag([100 100 -200]);

                % Compile physical metadata and the independently specified preparation
                spin_system=create(sys,inter);
                v=[0;1;1;0]/sqrt(2);
                observed=v*v';
                times=linspace(0,200e-6,31);

                % Independently construct the full dense Hamiltonian from the input Hz tensors
                dimensions=spin_system.comp.mults;
                dimension=prod(dimensions);
                operators=cell(numel(dimensions),3);
                h=zeros(dimension);
                for j=1:numel(dimensions)

                    % Local angular momentum matrices and their tensor-product embeddings
                    q=(dimensions(j)-1)/2;
                    m=q:-1:-q;
                    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
                    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
                    for a=1:3
                        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
                                            kron(local{a},eye(prod(dimensions(j+1:end)))));
                    end

                end
                for j=1:numel(dimensions)
                    for k=j:numel(dimensions)

                        % Spinach ordered pair matrices transpose into Cartesian operators
                        tensor=inter.coupling.matrix{j,k};
                        if isempty(tensor)
                            continue;
                        end
                        for a=1:3
                            for b=1:3
                                h=h+2*pi*tensor(b,a)*operators{j,a}*operators{k,b};
                            end
                        end

                    end
                end

                % Propagate an explicit product density and trace over the whole finite bath
                d_obs=size(observed,1);
                db=dimension/d_obs;
                rho0=kron(observed,eye(db)/db);
                reference=zeros(size(answer));
                for f=1:numel(times)

                    % This independent calculation uses expm, not the CCE eigensolver
                    u=expm(-1i*h*times(f));
                    rho=u*rho0*u';
                    reduced=zeros(d_obs);
                    for k=1:db
                        ids=k:db:dimension;
                        reduced=reduced+rho(ids,ids);
                    end
                    reference(f,:)=reduced(:).';

                end
                reference_error=max(abs(answer-reference),[],'all');

                % Preserve independent physical arrays for the stage-0 comparison
                snapshot=struct('a',a,'b',b,'db',db,'d_obs',d_obs,...
                                'dimension',dimension,'dimensions',dimensions,'h',h,...
                                'local',{local},'m',m,'operators',{operators},...
                                'plus',plus,'reduced',reduced,'reference',reference,...
                                'rho',rho,'rho0',rho0,'u',u,'v',v,'observed',observed,...
                                'reference_error',reference_error,'tensor',tensor,'ids',ids);

                % Retain this case reference in the public capture structure
                actual.reference_error=reference_error;
                reference_capture.gcce_highspin_siv=snapshot;
            case 'gcce_votpp'

                % Start fresh physical inputs for this independent reference
                sys=struct();
                inter=struct();

                % Magnet and isotopes; zero electron g isolates the model energy scales
                sys.magnet=0;
                sys.isotopes={'E','51V','1H'};
                sys.output='hush';
                sys.parallel={'local',1};
                inter.zeeman.matrix=repmat({zeros(3)},1,3);
                inter.coupling.matrix=cell(3);

                % Standard Spinach coupling and zero-field-splitting tensors in Hz
                inter.coupling.matrix{1,2}=diag([15 15 50])*1e6;
                inter.coupling.matrix{1,3}=diag([2000 3000 7000]);
                inter.coupling.matrix{2,3}=diag([100 100 -200]);
                inter.coupling.matrix{2,2}=diag([-1 -1 2])*1e5;

                % Compile physical metadata and the independently specified preparation
                spin_system=create(sys,inter);
                v=zeros(16,1);
                v([1 9])=1/sqrt(2);
                observed=v*v';
                times=linspace(0,200e-6,31);

                % Independently construct the full dense Hamiltonian from the input Hz tensors
                dimensions=spin_system.comp.mults;
                dimension=prod(dimensions);
                operators=cell(numel(dimensions),3);
                h=zeros(dimension);
                for j=1:numel(dimensions)

                    % Local angular momentum matrices and their tensor-product embeddings
                    q=(dimensions(j)-1)/2;
                    m=q:-1:-q;
                    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
                    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
                    for a=1:3
                        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
                                            kron(local{a},eye(prod(dimensions(j+1:end)))));
                    end

                end
                for j=1:numel(dimensions)
                    for k=j:numel(dimensions)

                        % Spinach ordered pair matrices transpose into Cartesian operators
                        tensor=inter.coupling.matrix{j,k};
                        if isempty(tensor)
                            continue;
                        end
                        for a=1:3
                            for b=1:3
                                h=h+2*pi*tensor(b,a)*operators{j,a}*operators{k,b};
                            end
                        end

                    end
                end

                % Propagate an explicit product density and trace over the whole finite bath
                d_obs=size(observed,1);
                db=dimension/d_obs;
                rho0=kron(observed,eye(db)/db);
                reference=zeros(size(answer));
                for f=1:numel(times)

                    % This independent calculation uses expm, not the CCE eigensolver
                    u=expm(-1i*h*times(f));
                    rho=u*rho0*u';
                    reduced=zeros(d_obs);
                    for k=1:db
                        ids=k:db:dimension;
                        reduced=reduced+rho(ids,ids);
                    end
                    reference(f,:)=reduced(:).';

                end
                reference_error=max(abs(answer-reference),[],'all');

                % Preserve independent physical arrays for the stage-0 comparison
                snapshot=struct('a',a,'b',b,'db',db,'d_obs',d_obs,...
                                'dimension',dimension,'dimensions',dimensions,'h',h,...
                                'local',{local},'m',m,'operators',{operators},...
                                'plus',plus,'reduced',reduced,'reference',reference,...
                                'rho',rho,'rho0',rho0,'u',u,'v',v,'observed',observed,...
                                'reference_error',reference_error,'tensor',tensor,'ids',ids);

                % Retain this case reference in the public capture structure
                actual.reference_error=reference_error;
                reference_capture.gcce_votpp=snapshot;
            case 'gcce_finite_drive'

                % Independent two-spin Hamiltonian including the finite rectangular drive
                x=[0 1;1 0]/2;
                y=[0 -1i;1i 0]/2;
                z=diag([1 -1])/2;
                v=[1;1]/sqrt(2);
                rho0=kron(v*v',eye(2)/2);
                h=2*pi*(1100*kron(x,x)+900*kron(y,y)+2200*kron(z,z));
                drive=4e4*kron(x,eye(2));
                times=diagnostics.parameters.times(:);
                exact=zeros(size(answer));
                for f=1:numel(times)

                    % Keep the independent expm calculation and its original arithmetic order
                    t=times(f);
                    duration=min(t,20e-6);
                    u=expm(-1i*h*(t-duration))*expm(-1i*(h+4e4*kron(x,eye(2)))*duration);
                    fullrho=u*rho0*u';
                    small=zeros(2);
                    for k=1:2
                        ids=k:2:4;
                        small=small+fullrho(ids,ids);
                    end
                    exact(f,:)=small(:).';

                end

                % Preserve the physical reference arrays retained by the stage-zero capture
                reference_error=max(abs(answer-exact),[],'all');
                populations=real(answer(:,1:3:end));
                snapshot=struct('drive',drive,'duration',duration,'exact',exact,...
                                'fullrho',fullrho,'h',h,'ids',ids,'populations',populations,...
                                'rho0',rho0,'small',small,'t',t,'u',u,'v',v,'x',x,'y',y,'z',z);

                % Retain this case reference in the public capture structure
                actual.reference_error=reference_error;
                reference_capture.gcce_finite_drive=snapshot;
            case 'gcce_two_centres'

                % Start fresh physical inputs for this independent reference
                sys=struct();
                inter=struct();
                ref_params=struct();

                % Magnet and isotopes; zero electron g isolates the model energy scales
                sys.magnet=0;
                sys.isotopes={'E','E','13C','13C'};
                sys.output='hush';
                sys.parallel={'local',1};
                inter.zeeman.matrix=repmat({zeros(3)},1,4);
                inter.coupling.matrix=cell(4);

                % Standard Spinach coupling and zero-field-splitting tensors in Hz
                inter.coupling.matrix{1,2}=diag([700 500 1200]);
                inter.coupling.matrix{1,3}=diag([1700 900 2000]);
                inter.coupling.matrix{2,4}=diag([1300 800 1600]);
                inter.coupling.matrix{3,4}=diag([100 100 -200]);

                % Compile metadata, without a full-bath basis
                spin_system=create(sys,inter);

                % Generalised CCE propagates the entire observed density matrix
                ref_params.observed_spins=[1 2];
                ref_params.bath_spins=[3 4];
                ref_params.method='gcce';
                ref_params.expansion='standard';
                ref_params.quantity='density';
                ref_params.order=2;
                v=[1;0;0;1]/sqrt(2);
                ref_params.rho0=v*v';
                ref_params.assumptions='labframe';
                ref_params.orientation=[0 0 0];
                ref_params.interaction_cutoff=Inf;
                ref_params.meanfield='none';
                ref_params.graph=struct('mode','complete');
                ref_params.times=linspace(0,200e-6,31);
                ref_params.bath_state=struct('mode','mixed_product');
                ref_params.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);

                % Run full order for the small bath, retaining all observed density-matrix off-diagonal elements
                ref_result=diagnostics.reduction;
                ref_result.value=answer;
                ref_result.times=ref_params.times(:);
                ref_result.density=reshape(answer.',sqrt(size(answer,2)),sqrt(size(answer,2)),[]);
                ref_result.diagnostics.trace_error=zeros(numel(ref_params.times),1);
                ref_result.diagnostics.minimum_eigenvalue=zeros(numel(ref_params.times),1);

                % Check physicality before evaluating spin correlations
                negativity=NaN(size(ref_result.times));
                for f=1:numel(ref_result.times)

                    % Partial transpose of the second observed electron, no eigenvalue clipping
                    rho=ref_result.density(:,:,f);
                    if all(ref_result.valid(f,:))&&abs(trace(rho)-1)<1e-8&&...
                       norm(rho-rho','fro')<1e-8&&min(eig((rho+rho')/2))>-1e-8
                        blocks=reshape(rho,2,2,2,2);
                        partial=reshape(permute(blocks,[3 2 1 4]),4,4);
                        negativity(f)=(sum(abs(eig(partial)))-1)/2;
                    end

                end
                ref_result.negativity=negativity;

                % Independently construct the full dense Hamiltonian from the input Hz tensors
                dimensions=spin_system.comp.mults;
                dimension=prod(dimensions);
                operators=cell(numel(dimensions),3);
                h=zeros(dimension);
                for j=1:numel(dimensions)

                    % Local angular momentum matrices and their tensor-product embeddings
                    q=(dimensions(j)-1)/2;
                    m=q:-1:-q;
                    plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
                    local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
                    for a=1:3
                        operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
                                            kron(local{a},eye(prod(dimensions(j+1:end)))));
                    end

                end
                for j=1:numel(dimensions)
                    for k=j:numel(dimensions)

                        % Spinach ordered pair matrices transpose into Cartesian operators
                        tensor=inter.coupling.matrix{j,k};
                        if isempty(tensor)
                            continue;
                        end
                        for a=1:3
                            for b=1:3
                                h=h+2*pi*tensor(b,a)*operators{j,a}*operators{k,b};
                            end
                        end

                    end
                end

                % Propagate an explicit product density and trace over the whole finite bath
                d_obs=size(ref_params.rho0,1);
                db=dimension/d_obs;
                rho0=kron(ref_params.rho0,eye(db)/db);
                reference=zeros(size(ref_result.value));
                for f=1:numel(ref_result.times)

                    % This independent calculation uses expm, not the CCE eigensolver
                    u=expm(-1i*h*ref_result.times(f));
                    rho=u*rho0*u';
                    reduced=zeros(d_obs);
                    for k=1:db
                        ids=k:db:dimension;
                        reduced=reduced+rho(ids,ids);
                    end
                    reference(f,:)=reduced(:).';

                end
                ref_result.reference_error=max(abs(ref_result.value-reference),[],'all');

                % Plot observed populations without silently repairing the density matrix
                populations=zeros(numel(ref_result.times),4);
                for f=1:numel(ref_result.times)

                    % Read the diagonal in the documented observed tensor-product order
                    populations(f,:)=real(diag(ref_result.density(:,:,f))).';

                end

                % Preserve every frozen physical reference array
                reference_error=ref_result.reference_error;
                snapshot=struct();
                snapshot.a=a;
                snapshot.b=b;
                snapshot.blocks=blocks;
                snapshot.d_obs=d_obs;
                snapshot.db=db;
                snapshot.dimension=dimension;
                snapshot.dimensions=dimensions;
                snapshot.h=h;
                snapshot.ids=ids;
                snapshot.local=local;
                snapshot.m=m;
                snapshot.negativity=negativity;
                snapshot.operators=operators;
                snapshot.partial=partial;
                snapshot.plus=plus;
                snapshot.populations=populations;
                snapshot.reduced=reduced;
                snapshot.reference=reference;
                snapshot.rho=rho;
                snapshot.rho0=rho0;
                snapshot.u=u;
                snapshot.v=v;

                % Retain this case reference in the public capture structure
                actual.reference_error=reference_error;
                reference_capture.gcce_two_centres=snapshot;
            case 'mecce_dissipative_bath'

                % Start fresh physical inputs for this independent reference
                ref_params=struct();

                % Full density ME-gCCE with an explicitly supplied nuclear jump rate
                ref_params.observed_spins=1;
                ref_params.bath_spins=2;
                ref_params.method='megcce';
                ref_params.expansion='standard';
                ref_params.quantity='density';
                ref_params.order=1;
                v=[1;1]/sqrt(2);
                ref_params.rho0=v*v';
                ref_params.assumptions='labframe';
                ref_params.orientation=[0 0 0];
                ref_params.interaction_cutoff=Inf;
                ref_params.meanfield='none';
                ref_params.graph=struct('mode','complete');
                ref_params.times=linspace(0,200e-6,31);
                ref_params.bath_state=struct('mode','mixed_product');
                ref_params.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',2^28);
                jump=[0 1;0 0];
                rate=3000;
                ref_params.jumps=struct('support',2,'operator',jump,'rate',rate);

                % Propagate the reduced density with the cluster implementation
                ref_result=diagnostics.reduction;
                ref_result.value=answer;
                ref_result.times=ref_params.times(:);
                ref_result.density=reshape(answer.',sqrt(size(answer,2)),sqrt(size(answer,2)),[]);
                ref_result.diagnostics.trace_error=zeros(numel(ref_params.times),1);
                ref_result.diagnostics.minimum_eigenvalue=zeros(numel(ref_params.times),1);

                % Independently construct the column-vectorised dense master equation
                x=[0 1;1 0]/2;
                y=[0 -1i;1i 0]/2;
                z=diag([1 -1])/2;
                h=2*pi*(1100*kron(x,x)+900*kron(y,y)+2200*kron(z,z));
                a=kron(eye(2),jump);
                q=a'*a;
                unit=eye(4);
                liouvillian=-1i*(kron(unit,h)-kron(h.',unit))+...
                    rate*(kron(conj(a),a)-(kron(unit,q)+kron(q.',unit))/2);
                rho0=kron(ref_params.rho0,eye(2)/2);
                reference=zeros(size(ref_result.value));
                for f=1:numel(ref_result.times)

                    % Dense exponential followed by an explicit nuclear partial trace
                    rho=reshape(expm(liouvillian*ref_result.times(f))*rho0(:),4,4);
                    reduced=zeros(2);
                    for k=1:2

                        % Tensor order is observed electron, then bath nucleus
                        ids=k:2:4;
                        reduced=reduced+rho(ids,ids);

                    end
                    reference(f,:)=reduced(:).';

                end
                ref_result.reference_error=max(abs(ref_result.value-reference),[],'all');

                % Preserve every frozen physical reference array
                reference_error=ref_result.reference_error;
                snapshot=struct();
                snapshot.a=a;
                snapshot.h=h;
                snapshot.ids=ids;
                snapshot.jump=jump;
                snapshot.liouvillian=liouvillian;
                snapshot.rate=rate;
                snapshot.reduced=reduced;
                snapshot.reference=reference;
                snapshot.rho=rho;
                snapshot.rho0=rho0;
                snapshot.unit=unit;
                snapshot.v=v;
                snapshot.x=x;
                snapshot.y=y;
                snapshot.z=z;

                % Retain this case reference in the public capture structure
                actual.reference_error=reference_error;
                reference_capture.mecce_dissipative_bath=snapshot;
                result=test_true(result,'Liouville arithmetic bounds',...
                                 all(diagnostics.reduction.abs_error_bound>0,'all'),...
                                 'The stock Liouville action retains its arithmetic estimate.');
        end
        reference_capture.(name).diagnostics=diagnostics;
    elseif strcmp(name,'cce_deer_bandwidth')
        [answer,diagnostics]=cases{n,2}();
        actual=diagnostics.reduction;
        actual.value=answer;
        actual.density=[];

        % Assemble the independent Cartesian Hamiltonian from physical toy couplings
        sx=[0 1;1 0]/2;
        sy=[0 -1i;1i 0]/2;
        sz=diag([0.5 -0.5]);
        x1=kron(sx,eye(4));
        x2=kron(eye(2),kron(sx,eye(2)));
        y2=kron(eye(2),kron(sy,eye(2)));
        z2=kron(eye(2),kron(sz,eye(2)));
        x3=kron(eye(4),sx);
        y3=kron(eye(4),sy);
        z3=kron(eye(4),sz);
        z1=kron(sz,eye(4));
        H=2*pi*(20000*z1*z2+12000*z1*z3+200*x2*x3+200*y2*y3-400*z2*z3);
        rho0=kron([0.5 0.5;0.5 0.5],eye(4)/4);
        detector=kron([0 1;0 0],eye(4));
        pump=expm(-1i*pi*x2);
        refocus=expm(-1i*pi*x1);
        reference=complex(zeros(size(answer)));
        times=diagnostics.parameters.times;

        % The pump precedes the refocusing pulse; zero time is the initial reference
        for f=1:numel(times)
            u=expm(-1i*H*times(f)/4);
            v=expm(-1i*H*times(f)/2);
            if times(f)==0
                rho=rho0;
            else
                propagator=v*refocus*u*pump*u;
                rho=propagator*rho0*propagator';
            end
            reference(f)=trace(detector'*rho)/trace(detector'*rho0);
        end
        reference_error=max(abs(answer-reference));
        snapshot=struct('reference',reference,'reference_error',reference_error);

        % Retain this case reference in the public capture structure
        actual.reference_error=reference_error;
        reference_capture.cce_deer_bandwidth=snapshot;
    elseif strcmp(name,'cce_molecular_powder')
        [answer,diagnostics]=cases{n,2}();
        actual=diagnostics.reduction;
        actual.value=answer;
        actual.density=[];

        % Start fresh physical inputs for this independent reference
        sys=struct();
        inter=struct();

        % Define the same physical tensors without a cluster-expansion service
        sys.magnet=0.01;
        sys.isotopes={'E','2H','2H'};
        sys.output='hush';
        sys.parallel={'local',1};
        inter.zeeman.matrix=repmat({zeros(3)},1,3);
        inter.coupling.matrix=cell(3);
        inter.coupling.matrix{1,2}=[0 0 0;0 0 0;40000 0 20000];
        inter.coupling.matrix{1,3}=[0 0 0;0 0 0;20000 0 -10000];
        inter.coupling.matrix{2,3}=diag([500 500 -1000]);
        inter.coupling.matrix{2,2}=diag([-1 -1 2])*2e5;
        inter.coupling.matrix{3,3}=diag([-1 -1 2])*5e4;
        spin_system=create(sys,inter);
        dimensions=spin_system.comp.mults;
        operators=cell(3,3);
        dimension=prod(dimensions);
        for j=1:3
            q=(dimensions(j)-1)/2;
            m=q:-1:-q;
            plus=diag(sqrt(q*(q+1)-m(2:end).*(m(2:end)+1)),1);
            local={(plus+plus')/2,(plus-plus')/(2i),diag(m)};
            for a=1:3
                operators{j,a}=kron(eye(prod(dimensions(1:j-1))),...
                                    kron(local{a},eye(prod(dimensions(j+1:end)))));
            end
        end

        % Rotate every tensor and contract it with explicit Cartesian spin matrices
        grid=diagnostics.grid;
        times=diagnostics.parameters.times;
        values=complex(zeros(numel(times),numel(grid.weights)));
        rho0=kron(ones(2)/2,eye(dimension/2)/(dimension/2));
        pulse=expm(-1i*pi*operators{1,1});
        for f=1:numel(grid.weights)
            rotation=euler2dcm(grid.alphas(f),grid.betas(f),grid.gammas(f));
            h=zeros(dimension);
            for j=1:3
                tensor=rotation*spin_system.inter.zeeman.matrix{j}*rotation';
                for a=1:3
                    h=h+tensor(3,a)*operators{j,a};
                end
                for k=j:3
                    tensor=spin_system.inter.coupling.matrix{j,k};
                    if isempty(tensor)
                        continue;
                    end

                    % Ordered Spinach pair cells contract in transposed Cartesian order
                    if k~=j
                        tensor=tensor.';
                    end
                    tensor=rotation*tensor*rotation';
                    for a=1:3
                        for b=1:3
                            h=h+tensor(a,b)*operators{j,a}*operators{k,b};
                        end
                    end
                end
            end

            % Full-bath exponentials do not use the CCE eigensolver or its reduction
            for k=1:numel(times)
                u=eye(dimension);
                if times(k)>0
                    free=expm(-1i*h*times(k)/2);
                    u=free*pulse*free;
                end
                rho=u*rho0*u';
                observed=zeros(2);
                for j=1:dimension/2
                    indices=j:dimension/2:dimension;
                    observed=observed+rho(indices,indices);
                end
                values(k,f)=2*observed(1,2);
            end
        end
        reference=values*grid.weights;
        reference_error=max(abs(reference-answer));
        snapshot=struct('reference',reference,'orientation_values',values,'reference_error',reference_error);

        % Retain this case reference in the public capture structure
        actual.reference_error=reference_error;
        reference_capture.cce_molecular_powder=snapshot;
    elseif strcmp(name,'cce_polarized_bath')
        [answer,diagnostics]=cases{n,2}();
        actual=diagnostics.reduction;
        actual.value=answer;

        % Start fresh physical inputs for this independent reference
        sys=struct();
        inter=struct();

        % Standard Spinach magnet, isotope list, and illustrative Hz tensors
        sys.magnet=0.01;
        sys.isotopes={'E','1H','1H'};
        sys.output='hush';
        sys.parallel={'local',1};
        inter.zeeman.matrix={zeros(3),zeros(3),zeros(3)};
        inter.coupling.matrix=cell(3);
        inter.coupling.matrix{1,2}=[0 0 0;0 0 0;15000 0 8000].';
        inter.coupling.matrix{1,3}=[0 0 0;0 0 0;9000 0 12000].';
        inter.coupling.matrix{2,3}=diag([400 400 -800]);
        spin_system=create(sys,inter);

        % The same physical populations, averaged exhaustively rather than sampled
        ref_params=diagnostics.parameters;
        population=[0.8;0.2];
        ref_params.cce.bath_state=struct('mode','explicit_product',...
                                         'rho',{{diag(population),diag(population)}});
        [exact,exact_diagnostics]=clexpn(spin_system,@hahn_echo,ref_params,'labframe');
        sampling_error=max(abs(answer-exact));
        snapshot=struct('population',population,'exact',exact,...
                        'exact_diagnostics',exact_diagnostics,'sampling_error',sampling_error);

        % Retain this case reference in the public capture structure
        actual.sampling_error=sampling_error;
        reference_capture.cce_polarized_bath=snapshot;
        actual.density=[];
        reference_capture.cce_polarized_bath.diagnostics=diagnostics;
        parts=diagnostics.reduction.samples;
        weights=diagnostics.reduction.weights;
        bound=zeros(size(answer));
        magnitude=zeros(size(answer));
        for sample=1:numel(parts)
            bound=bound+weights(sample)*parts{sample}.abs_error_bound;
            magnitude=magnitude+weights(sample)*abs(parts{sample}.value);
        end
        operations=4*numel(parts);
        rounding=operations*eps/(1-operations*eps);
        result=test_close(result,'Sampled bath arithmetic bound',...
                          diagnostics.reduction.abs_error_bound,bound+rounding*magnitude,0,0,...
                          'Deterministic summation error is distinct from estimator uncertainty.');
    else
        actual=cases{n,2}();
    end
    result=test_true(result,[name ' valid'],all(actual.valid,'all'),...
                     'Every displayed time point has a valid expansion.');
    result=test_true(result,[name ' finite'],all(isfinite(actual.value),'all'),...
                     'Every returned complex value is finite.');
    if isfield(actual,'reference_error')

        % Independent small-system targets use an absolute trace-error criterion
        result=test_close(result,[name ' reference'],actual.reference_error,...
                          0,1e-7,0,'Analytic, dense expm, or pinned Python reference.');

    end
    if ~isempty(actual.density)

        % Full-order finite baths must preserve observed density physicality
        result=test_true(result,[name ' trace'],...
                         max(actual.diagnostics.trace_error)<1e-8,...
                         'Reduced density matrices retain unit trace.');
        result=test_true(result,[name ' positive'],...
                         min(actual.diagnostics.minimum_eigenvalue)>-1e-8,...
                         'No material negative population at full bath order.');

    end
    if strcmp(name,'cce_deer_bandwidth')

        % Excluding the pump resonance must remove its selected pulse
        parameters=diagnostics.parameters;
        parameters.pump.center=1e7;
        sys.magnet=0.01;
        sys.isotopes={'E','E','E'};
        sys.output='hush';
        sys.parallel={'local',1};
        inter=struct('zeeman',struct('scalar',{{0,0,0}}));
        spin_system=create(sys,inter);
        unpumped=clexpn(spin_system,@bath_pumped_hahn,parameters,'labframe');
        parameters.pump.allow_empty=false;
        rejected=false;
        try
            clexpn(spin_system,@bath_pumped_hahn,parameters,'labframe');
        catch exception
            rejected=strcmp(exception.identifier,'bath_pumped_hahn:empty');
        end
        result=test_true(result,'DEER example excludes pump resonance',...
                         rejected&&max(abs(unpumped-1))<1e-12,...
                         'The native selector removes the pulse and enforces allow_empty.');

    end
    if strcmp(name,'cce_molecular_powder')

        % Completed orientation expansions must reproduce the powder average
        values=cat(2,diagnostics.values{:});
        result=test_close(result,'Powder weighted response',actual.value,...
                          values*diagnostics.grid.weights,2e-15,2e-15,...
                          'Powder averages completed physical responses.');
        result=test_close(result,'Powder final order response',...
                          actual.by_order(:,:,end),actual.value,0,0,...
                          'The final expansion order equals the returned signal.');
        result=test_true(result,'Powder finite diagnostics',...
                         isequal(size(actual.abs_error_bound),size(actual.value))&&...
                         all(isfinite(actual.abs_error_bound)),...
                         'Numerical bounds retain the powder output shape.');
        expected=zeros(size(actual.abs_error_bound));
        valid=true(size(actual.valid));
        for k=1:numel(diagnostics.grid.weights)
            current=diagnostics.orientations{k}.reduction;
            expected=expected+diagnostics.grid.weights(k)*current.abs_error_bound;
            valid=valid&current.valid;
        end
        result=test_close(result,'Powder averaged numerical bounds',...
                          actual.abs_error_bound,expected,0,0,...
                          'Positive quadrature weights also average absolute bounds.');
        result=test_true(result,'Powder combined numerical validity',...
                         isequal(actual.valid,valid),...
                         'Every contributing orientation must satisfy numerical validity.');
        result=test_close(result,'Powder zero-time analytic limit',...
                          actual.value(1),1,2e-15,0,...
                          'The normalized initial Hahn signal is unity.');

    end
    fprintf('Example %s completed in %.3f s.\n',name,toc(started));
    close all;

end

end


