% ME-CCE and ME-gCCE against closed and independent dense Liouville limits.
%
% ilya.kuprov@weizmann.ac.il

function result=test_mecce_exact()

% Construct a central spin and two interacting bath spins
result=new_test_result('kernel/mecce_exact','ME exact dynamics',...
                       'Zero jumps, dense Liouville full bath and correlated support.');
sys.isotopes={'1H','1H','1H'}; sys.magnet=0; sys.output='hush';
inter.zeeman.scalar={0,0,0}; inter.coupling.scalar=cell(3);
s=create(sys,inter); s.inter.zeeman.matrix{2}(3,3)=3;
s.inter.coupling.matrix{1,2}=diag([0.7 0.4 1.1]);
s.inter.coupling.matrix{2,3}=diag([0.3 0.5 0.9]);
p=parameters(); p.times=linspace(0,0.8,15);
p.central=struct('alpha',[1;0],'beta',[0;1],'rho',[1 1;1 1]/2);
p.jumps=struct('support',{},'operator',{},'rate',{});

% Compare selective explicit-generator ME pulses with closed and dense evolution
q=p; spin=pauli(2); q.method='gcce'; q.quantity='density'; q.times=[0 0.07 0.2];
q.sequence=struct('type','ideal','target',2,'placement','fraction',...
                  'at',0.4,'operator',spin.x,'angle',pi);
closed_pulse=cce(s,q); q.method='megcce'; opened_pulse=cce(s,q);
result=test_close(result,'Selective operator ME pulse matches closed evolution',...
                  opened_pulse.density,closed_pulse.density,2e-11,2e-11,...
                  'Zero-jump master evolution must use the same retained-site dimensions as the Hamiltonian.');
q.method='gcce'; model=cce_prepare(s,q); H=cce_hamiltonian(model,[2 3],struct());
spin=pauli(2); pulse=kron(eye(2),kron(expm(-1i*pi*spin.x),eye(2)));
initial=kron(q.central.rho,eye(4)/4); pulse_error=0;
for k=1:numel(q.times)
    t=q.times(k); U=expm(-1i*H*0.6*t)*pulse*expm(-1i*H*0.4*t);
    if t==0, U=eye(8); end
    full_density=U*initial*U'; reduced=zeros(2);
    for a=1:2
        for b=1:2
            reduced(a,b)=trace(full_density((a-1)*4+(1:4),(b-1)*4+(1:4)));
        end
    end
    result=test_close(result,'Selective operator ME pulse dense reference',...
                      opened_pulse.density(:,:,k),reduced,2e-11,2e-11,...
                      'An independent Kronecker pulse embeds the bath generator at physical site two.');
    pulse_error=max(pulse_error,norm(opened_pulse.density(:,:,k)-reduced,'fro'));
end
fprintf('Selective ME pulse dense discrepancy %.3g.\n',pulse_error);

% Require zero-jump trajectories to retain complex closed-system phases
p.method='cce'; closed=cce(s,p); p.method='mecce'; opened=cce(s,p);
result=test_close(result,'C7 numerical check',opened.value,closed.value,2e-11,2e-11,'Conditional zero jumps.');
p.method='gcce'; p.quantity='density'; closed=cce(s,p);
p.method='megcce'; opened=cce(s,p);
result=test_close(result,'C7 numerical check',opened.density,closed.density,2e-11,2e-11,'Generalized zero jumps.');

% Preserve raw coherence unless normalization was explicitly requested
q=p; q.quantity='coherence'; q.central.alpha=[1;0]; q.central.beta=[0;1];
q.method='gcce'; closed_coherence=cce(s,q); q.method='megcce'; raw_coherence=cce(s,q);
result=test_close(result,'ME raw coherence default',raw_coherence.value,closed_coherence.value,...
                  2e-11,2e-11,'An empty jump array must not normalize raw generalized coherence.');
q.normalized=false; explicit_raw=cce(s,q);
result=test_close(result,'ME explicit raw coherence',explicit_raw.value,closed_coherence.value,...
                  2e-11,2e-11,'Explicit false and the omitted normalization option must agree.');
q.normalized=true; normalized_coherence=cce(s,q);
result=test_close(result,'ME explicit normalized coherence',normalized_coherence.value,...
                  2*closed_coherence.value,2e-11,2e-11,...
                  'Explicit normalization divides by the balanced initial coherence of one half.');
fprintf('R1230 closed %.16g; ME default %.16g.\n',...
        real(closed_coherence.value(1)),real(raw_coherence.value(1)));

% Refuse a multi-site operator instead of silently applying identity
q=p; q.method='megcce'; q.quantity='density'; q.central.rho=diag([1 0]);
q.sequence=struct('type','ideal','placement','fraction','at',0.5,'target',[1 2],...
                  'operator',kron([0 1;1 0]/2,eye(2)),'angle',pi);
rejected=false;
try
    cce(s,q);
catch err
    rejected=strcmp(err.identifier,'cce_master:operator_target');
end
result=test_true(result,'ME multi-site operator rejection',rejected,...
                 'Unsupported joint central-bath generators must fail explicitly, never act as identity.');

% Validate whole-sequence Hahn timing and a zero-duration initial state
p.sequence=struct('type','ideal','placement','fraction','at',0.5,...
                  'target',1,'axis','x','angle',pi);
p.method='gcce'; closed_echo=cce(s,p); p.method='megcce'; opened_echo=cce(s,p);
result=test_close(result,'ME-gCCE Hahn zero jumps',opened_echo.density,closed_echo.density,...
                  2e-11,2e-11,'Physical pulse, whole-sequence clock and zero-time identity.');
p.method='cce'; p.quantity='coherence'; closed_echo=cce(s,p);
p.method='mecce'; opened_echo=cce(s,p);
result=test_close(result,'ME-CCE Hahn zero jumps',opened_echo.value,closed_echo.value,...
                  2e-11,2e-11,'Conditional branch swapping retains the complex trajectory.');
p.method='megcce'; p.quantity='density'; p.sequence=struct([]);

% Build an independent Pauli-product Liouville reference with bath relaxation
x=[0 1;1 0]/2; y=[0 -1i;1i 0]/2; z=diag([1 -1])/2; e=eye(2);
H=3*kron(e,kron(z,e))+0.7*kron(x,kron(x,e))+0.4*kron(y,kron(y,e))+...
  1.1*kron(z,kron(z,e))+0.3*kron(e,kron(x,x))+...
  0.5*kron(e,kron(y,y))+0.9*kron(e,kron(z,z));
j=[0 1;0 0]; L=kron(e,kron(j,e)); rate=0.6;
G=-1i*(kron(eye(8),H)-kron(H.',eye(8)))+...
  rate*(kron(conj(L),L)-(kron(eye(8),L'*L)+kron((L'*L).',eye(8)))/2);
p.jumps=struct('support',2,'operator',j,'rate',rate); opened=cce(s,p);
rho=kron(p.central.rho,eye(4)/4); reference=zeros(2,2,numel(p.times));
for k=1:numel(p.times)
    evolved=reshape(expm(G*p.times(k))*rho(:),8,8);
    for a=1:2
        for b=1:2
            reference(a,b,k)=trace(evolved((a-1)*4+(1:4),(b-1)*4+(1:4)));
        end
    end
end
result=test_close(result,'C7 numerical check',opened.density,reference,2e-11,2e-11,'Independent full Liouville density.');

% Validate the full conditional bath block against an independent Liouvillian
ha=H(1:4,1:4); hb=H(5:8,5:8); lb=kron(j,eye(2));
generator=-1i*(kron(eye(4),ha)-kron(hb.',eye(4)))+...
          rate*(kron(conj(lb),lb)-(kron(eye(4),lb'*lb)+kron((lb'*lb).',eye(4)))/2);
p.method='mecce'; p.quantity='coherence'; conditional=cce(s,p); rb=eye(4)/4;
conditional_reference=zeros(numel(p.times),1);
for k=1:numel(p.times)
    conditional_reference(k)=trace(reshape(expm(generator*p.times(k))*rb(:),4,4));
end
result=test_close(result,'Full-bath conditional ME reference',conditional.value,conditional_reference,...
                  2e-11,2e-11,'Independent column-vectorized left/right Lindblad blocks.');
p.method='megcce'; p.quantity='density';

% Compare a dissipative echo directly with independently exponentiated blocks
p.sequence=struct('type','ideal','placement','fraction','at',0.5,...
                  'target',1,'axis','x','angle',pi);
echo=cce(s,p); pulse=kron(-2i*x,eye(4));
for k=2:numel(p.times)
    half=expm(G*p.times(k)/2); density=reshape(half*rho(:),8,8);
    density=pulse*density*pulse'; density=reshape(half*density(:),8,8);
    reduced=zeros(2);
    for a=1:2
        for b=1:2
            reduced(a,b)=trace(density((a-1)*4+(1:4),(b-1)*4+(1:4)));
        end
    end
    result=test_close(result,'Dissipative Hahn reference',echo.density(:,:,k),reduced,...
                      2e-11,2e-11,'Independent dense Liouville evolution on both sides of the pulse.');
end
p.sequence=struct([]);

% Refuse a split correlated jump and accept its complete ordered support
p.jumps=struct('support',[3 2],'operator',kron(j,z),'rate',rate);
m=cce_prepare(s,p); refused=false;
try
    cce_liouvillian(m,2,struct(),p);
catch err
    refused=strcmp(err.identifier,'cce:jumpSupport');
end
result=test_true(result,'C7 contract check',refused,'Split correlated jump is refused.');
[G,~]=cce_liouvillian(m,[2 3],struct(),p); L=kron(e,kron(z,j));
independent=-1i*(kron(eye(8),H)-kron(H.',eye(8)))+...
            rate*(kron(conj(L),L)-(kron(eye(8),L'*L)+kron((L'*L).',eye(8)))/2);
result=test_close(result,'C7 numerical check',G,independent,2e-12,2e-12,'Reversed support tensor ordering.');

% Compare the reconstructed PyCCE six-electron ME tutorial with declared readout
result=pycce_tutorial(result);

% Report numerical metrics through the standard harness without extra fields
fprintf('ME full-density max error %.3g; correlated generator %.3g\n',...
        max(abs(opened.density-reference),[],'all'),norm(G-independent,'fro'));

end

% Declare all physical and numerical choices explicitly
function p=parameters()
p=struct('central_spins',1,'bath_spins',[2 3],'method','mecce','quantity','coherence',...
         'expansion','standard','order',2,'times',0,'assumptions','labframe',...
         'orientation',[0 0 0],'interaction_cutoff',Inf,'meanfield','none');
p.graph=struct('mode','complete'); p.bath_state=struct('mode','mixed_product');
p.sequence=struct([]);
p.numerics=struct('atol',1e-8,'rtol',1e-8,'max_memory_bytes',1e9);
end

% Reconstructed PyCCE d74d3cb683bf3caa07231654314e554f440833b4: seed 2, six electrons, 300 G, not inferred from NPZ
function result=pycce_tutorial(result)

% Rebuild the physical Hamiltonian from explicit rad/s tensors
sys.isotopes=repmat({'1H'},1,7); sys.magnet=0; sys.output='hush';
inter.zeeman.scalar=num2cell(zeros(1,7)); inter.coupling.scalar=cell(7);
s=create(sys,inter); hyperfine=cell(1,6);
hyperfine{1}=[3.9420192177386513 -22.116700630680850 34.487995976573679;-22.116700630680850 11.411183737505917 29.152644214055385;34.487995976573679 29.152644214055385 -15.353202955244585];
hyperfine{2}=[-28.292250697131941 12.803178330108064 29.233624456667787;12.803178330108064 20.762467763179441 -7.1708294592669652;29.233624456667787 -7.1708294592669652 7.5297829339524966];
hyperfine{3}=[1.9964913747647539 -15.479148588312755 -7.8301813701485887;-15.479148588312755 -9.2070544257849889 -11.160828486911058;-7.8301813701485887 -11.160828486911058 7.2105630510202365];
hyperfine{4}=[22.619007901075260 -8.6653596795542036 48.699977138307304;-8.6653596795542036 41.101375134206762 19.261121031726550;48.699977138307304 19.261121031726550 -63.720383035282005];
hyperfine{5}=[220.17856332914070 294.90036908815011 -134.50201830526140;294.90036908815011 -403.54840800974375 338.00052908324028;-134.50201830526140 338.00052908324028 183.36984468060282];
hyperfine{6}=[14.988123972553279 -12.927663731394146 -5.0629406044633933;-12.927663731394146 -26.700547564036675 -17.769351941854591;-5.0629406044633933 -17.769351941854591 11.712423591483407];
for k=1:6
    s.inter.zeeman.matrix{k+1}(3,3)=5282579115.0000000;
    s.inter.coupling.matrix{1,k+1}=hyperfine{k};
end

% Retain every independently regenerated bath dipolar pair
s.inter.coupling.matrix{2,3}=[18.490820824481172 69.736420167682269 -19.807027525309429;69.736420167682269 -66.218949829561666 35.204003505681300;-19.807027525309429 35.204003505681300 47.728129005080490];
s.inter.coupling.matrix{2,4}=[11.777848557501965 -1.8687578090504040 -4.1412051266677494;-1.8687578090504040 6.1745250133336418 -13.671483566690702;-4.1412051266677494 -13.671483566690702 -17.952373570835604];
s.inter.coupling.matrix{2,5}=[1.4188198395876055 -2.2519617083688344 5.6702623230819071;-2.2519617083688344 3.2974163741645843 3.7786554276401856;5.6702623230819071 3.7786554276401856 -4.7162362137521896];
s.inter.coupling.matrix{2,6}=[10.489522578329916 -8.5931380214102262 9.0713357153134684;-8.5931380214102262 -4.1790738626647030 19.668655961366063;9.0713357153134684 19.668655961366063 -6.3104487156652196];
s.inter.coupling.matrix{2,7}=[1.7486674704754059 -6.9153803136068666 1.3559080209746208;-6.9153803136068666 -6.8599333630258563 2.4410531075557125;1.3559080209746208 2.4410531075557125 5.1112658925504517];
s.inter.coupling.matrix{3,4}=[6.7560320560754690 2.1891127602816001 1.8221031107458450;2.1891127602816001 -5.2870391595939354 -10.344958444830219;1.8221031107458450 -10.344958444830219 -1.4689928964815311];
s.inter.coupling.matrix{3,5}=[-2.2788583476396278 0.48192401682408326 7.3686590462854920;0.48192401682408326 4.8922097157370921 -0.49298636648361910;7.3686590462854920 -0.49298636648361910 -2.6133513680974634];
s.inter.coupling.matrix{3,6}=[-16.459727638010278 -8.6953115099296312 38.423992890366939;-8.6953115099296312 24.588765951888167 7.8036408132509383;38.423992890366939 7.8036408132509383 -8.1290383138778939];
s.inter.coupling.matrix{3,7}=[-9.6870055814039944 -12.870637455986344 1.5940709122087282;-12.870637455986344 0.54312498463924597 1.0818267917313946;1.5940709122087282 1.0818267917313946 9.1438805967647525];
s.inter.coupling.matrix{4,5}=[-2.5744373240231266 -10.936651578652313 2.6937767806548014;-10.936651578652313 -4.2605533646739469 2.9094197049211914;2.6937767806548014 2.9094197049211914 6.8349906886970730];
s.inter.coupling.matrix{4,6}=[5.1197832662530116 -8.2670172525041927 -1.9897261972521421;-8.2670172525041927 -12.332064500854594 -4.9932296429744696;-1.9897261972521421 -4.9932296429744696 7.2122812346015763];
s.inter.coupling.matrix{4,7}=[1.0048036705024974 -2.0332371091237604 -0.90953034204584016;-2.0332371091237604 -2.1684481765543482 -1.8634311542364834;-0.90953034204584016 -1.8634311542364834 1.1636445060518499];
s.inter.coupling.matrix{5,6}=[-5.8046301531264417 27.495265730264840 70.059427541647182;27.495265730264840 38.853722405872212 -33.356320933098807;70.059427541647182 -33.356320933098807 -33.049092252745766];
s.inter.coupling.matrix{5,7}=[8.6649131660432825 0.57124975950233725 0.75148052899144857;0.57124975950233725 -0.84591497551107175 -12.556493477788264;0.75148052899144857 -12.556493477788264 -7.8189981905322101];
s.inter.coupling.matrix{6,7}=[11.591666809583273 -27.291294263188128 -23.825767121862537;-27.291294263188128 -10.620517270242312 -35.418869760150514;-23.825767121862537 -35.418869760150514 -0.97114953934097914];

% Use the native rates and the exact whole-echo tutorial times
p=parameters(); p.bath_spins=2:7; p.times=linspace(0,0.001,11);
p.central=struct('alpha',[1;0],'beta',[0;1],'rho',[1 1;1 1]/2);
p.sequence=struct('type','ideal','placement','fraction','at',0.5,...
                  'target',1,'axis','x','angle',pi);
p.numerics=struct('atol',1e-3,'rtol',1e-3,'max_memory_bytes',1e9,'propagator','stock');
p.jumps=struct('support',{},'operator',{},'rate',{});
for k=2:7
    p.jumps(end+1)=struct('support',k,'operator',[0 1;0 0],'rate',500);
    p.jumps(end+1)=struct('support',k,'operator',[0 0;1 0],'rate',500);
end
r=cce(s,p);
reference=[1.0000000000000004;0.99999923394307144;0.99999409440779452;0.99998078627671738;0.99995608154472349;0.99991725565690237;0.99986203036107879;0.99978852242935690;0.99969519766699644;0.99958082968503426;0.99944446296394518]+...
          1i*[8.0506813103672611e-17;-1.4845050358724403e-11;-3.0705481891667437e-11;-4.7507112746376136e-11;-6.5180890168568098e-11;-8.3662944693980344e-11;-1.0289507687643922e-10;-1.2282225090608725e-10;-1.4339388638076503e-10;-1.6456348260287508e-10;-1.8628754486596782e-10];
result=test_close(result,'PyCCE reconstructed ME-CCE2',r.value,reference,5e-8,0,...
                  'Regenerated seeded coordinates, tensors, rates and Hahn sequence.');
conditional_error=max(abs(r.value-reference));

% Match PyCCE run/gcce.py:116-131 pulse-rotated final readout states
s.inter.zeeman.matrix{1}(3,3)=5282579115.0000000;
p.method='megcce'; p.central.alpha=[0;1]; p.central.beta=[1;0];
p.normalized=true; r=cce(s,p);
reference=[1.0000000000000020;0.99999696702632823;0.99996577952531562;0.99988384475122793;0.99978528981565296;0.99974208884696936;0.99976927314393593;0.99975844441771289;0.99954179272729826;0.99906118745487882;0.99848426483783148]+...
          1i*[-1.5501238790078434e-17;-0.0000045301957227707227;-0.000022166638316224802;-0.000022200000381094418;0.000036116963262306526;0.00012774472614107248;0.00014770063943958262;0.0000075750698896664989;-0.00023796879104827002;-0.00038050310796545583;-0.00022067377857413301];
result=test_close(result,'PyCCE reconstructed ME-gCCE2',r.value,reference,5e-8,0,...
                  'Final alpha/beta exchanged explicitly to match the PyCCE Hahn readout.');
result=test_true(result,'ME tutorial validity',all(r.valid),'Conservative lab-frame arithmetic estimates remain within declared budgets.');
fprintf('PyCCE ME-CCE2 %.3g; ME-gCCE2 %.3g; max bound %.3g\n',...
        conditional_error,max(abs(r.value-reference)),max(r.abs_error_bound));

end


