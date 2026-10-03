% Tests Spinach-facing JCAMP export of native magnetic resonance results.
% Syntax:
%
%                    result=test_jcamp_spinach()
%
% Outputs:
%
%    result - checks for dimensions, quadratures, coordinates, technique
%             metadata, scanner results, and refusal of ambiguous shapes
%
% talos@spindynamics.org

function result=test_jcamp_spinach()

% Check consistency
grumble();

% Prepare the conventional Spinach parameter and metadata structures
result=new_test_result('interfaces/jcamp_spinach','Spinach JCAMP wrappers',...
                       'native array layouts, quadratures, and physical coordinates survive export.');
spin_system.inter.magnet=14.1;
parameters=struct('spins',{{'1H'}},'sweep',800,'offset',100,'npoints',5,'zerofill',5);
info=struct('title','Pulse acquire','origin','Spinach','owner','Test',...
            'filename','','sequence','Pulse acquire','metadata',{cell(0,2)},...
            'delay',[0 0],'acquisition','SIMULTANEOUS');
fid=[1+2i;-3+4i;5-6i;-7-8i;9+10i];
text=jcamp_nmr(spin_system,parameters,fid,{'time'},info);
pairs=read_pairs(text);
result=test_true(result,'1D complex FID',isequal(pairs,[(0:4)'/800 real(fid); (0:4)'/800 imag(fid)]),...
                 'dwell, real channel, and imaginary channel retain their original values');
observe=abs(spin('1H')*14.1/(2*pi))*1e-6;
result=test_true(result,'observation frequency',contains(text,['##.OBSERVE FREQUENCY=' sprintf('%.17g',observe)]),...
                 'frequency is derived from the isotope and field in MHz');

% Use odd and even FFT grids without adding a duplicated edge frequency
for n=[5 6]
    spectrum=(1:n)'; text=jcamp_nmr(spin_system,parameters,spectrum,{'frequency'},info);
    expected=100+(-floor(n/2):ceil(n/2)-1)'*800/n;
    result=test_true(result,['FFT grid ' num2str(n)],isequal(read_pairs(text),[expected spectrum]),...
                     'frequency coordinates match the shifted FFT bins and Spinach plotting convention');
end

% Preserve the [F2,F1] shape and both States components of a 2D FID
parameters.spins={'13C','1H'}; parameters.sweep=[200 800]; parameters.offset=[20 100];
values=reshape(1:12,4,3)+1i*reshape(13:24,4,3);
fid=struct('pos',values,'neg',-values);
text=jcamp_nmr(spin_system,parameters,fid,{'time','time'},info);
pairs=read_pairs(text); expected=[];
for channel=[1 -1]
    for n=1:3
        expected=[expected; (0:3)'/800 channel*real(values(:,n));...
                  (0:3)'/800 channel*imag(values(:,n))]; %#ok<AGROW>
    end
end
result=test_true(result,'2D States amplitudes',isequal(pairs,expected),...
                 'pos and neg are separate LINK blocks without conjugation or recombination');
result=test_true(result,'2D dimension mapping',contains(text,'##VAR_NAME=F2 1H time, F1 13C time, Real signal, Imaginary signal')&&...
                 contains(text,'##PAGE=X2=0.005')&&contains(text,'##$SPINACH COMPONENT=neg'),...
                 'rows map to F2 and columns to F1, with component identity retained');
fid=struct('cos',values,'sin',2*values);
text=jcamp_nmr(spin_system,parameters,fid,{'time','time'},info);
result=test_true(result,'States cosine/sine',contains(text,'##TITLE=Pulse acquire/cos')&&...
                 contains(text,'##TITLE=Pulse acquire/sin'),...
                 'NOESY cosine and sine channels remain separately named');

% Accept the scalar homonuclear specifications used by COSY examples
parameters.spins={'1H'}; parameters.sweep=800; parameters.offset=100;
text=jcamp_nmr(spin_system,parameters,real(values),{'frequency','frequency'},info);
pairs=read_pairs(text);
result=test_true(result,'homonuclear 2D grid',isequal(pairs(:,1),repmat(ft_axis(100,800,4).',3,1)),...
                 'scalar nuclear, offset, and sweep specifications apply to both physical dimensions');

% Preserve mixed domains without mislabelling a tabulated frequency as time
text=jcamp_nmr(spin_system,parameters,values,{'time','frequency'},info);
result=test_true(result,'mixed-domain NMR',contains(text,'##DATA TYPE=NMR SPECTRUM')&&...
                 contains(text,'##UNITS=HZ, SECONDS')&&contains(text,'##.ACQUISITION MODE=SIMULTANEOUS')&&...
                 contains(text,'##.DELAY=(0, 0)'),...
                 'tabulated-axis type and fixed time coordinates describe a partly transformed array');

% Preserve the [F1,F2,F3] shape used by protein triple-resonance examples
parameters.spins={'15N','13C','1H'}; parameters.sweep=[100 200 800]; parameters.offset=[10 20 100];
values=reshape(1:24,2,3,4)+1i*reshape(25:48,2,3,4);
fid=struct('pos_pos',values,'pos_neg',2*values,'neg_pos',3*values,'neg_neg',4*values);
text=jcamp_nmr(spin_system,parameters,fid,{'time','time','time'},info);
pairs=read_pairs(text); expected=[];
for channel=1:4
    for n=1:12
        expected=[expected; (0:1)'/100 channel*real(values((n-1)*2+(1:2))).';...
                  (0:1)'/100 channel*imag(values((n-1)*2+(1:2))).']; %#ok<AGROW>
    end
end
result=test_true(result,'3D four-component amplitudes',isequal(pairs,expected),...
                 'all four HNCO/HNCA branches and every tensor sample retain their original values');
result=test_true(result,'3D coordinates',contains(text,['##PAGE=X2=' sprintf('%.17g',0.005) ', X3=' sprintf('%.17g',0.00125)])&&...
                 contains(text,'##$AXIS NUCLEI=15N, 13C, 1H'),...
                 'F2 and F3 coordinates follow MATLAB tensor indices and physical isotope order');
values=reshape(1:60,3,4,5);
text=jcamp_nmr(spin_system,parameters,values,{'frequency','frequency','frequency'},info);
result=test_true(result,'3D spectrum',contains(text,'##DATA TYPE=NMR SPECTRUM')&&...
                 isequal(read_pairs(text),[repmat(ft_axis(10,100,3).',20,1) values(:)]),...
                 'a processed triple-resonance spectrum follows the plot_3d layout');

% Export native field-swept and ENDOR row outputs without axis fitting
emr=struct('title','EMR result','origin','Spinach','owner','Test','filename','',...
           'detection','CW','method','SPECTRUM','description','Test physical parameters',...
           'metadata',{cell(0,2)});
parameters=struct('b_axis',[0.3 0.31 0.333 0.34],'mw_freq',9e9);
scan=[1+2i -3+4i 5-6i -7-8i];
text=jcamp_epr(spin_system,parameters,scan,'field',emr);
result=test_true(result,'field-swept EPR',isequal(read_pairs(text),...
                 [parameters.b_axis.' real(scan).'; parameters.b_axis.' imag(scan).'])&&...
                 contains(text,'##.MICROWAVE FREQUENCY=9000000000'),...
                 'returned irregular tesla coordinates and complex amplitudes are retained');
parameters.n_frq=[-20e6 -2e6 3e6 40e6]; emr.detection='PULSE'; emr.method='ENDOR';
text=jcamp_epr(spin_system,parameters,real(scan),'endor',emr);
result=test_true(result,'ENDOR scan',isequal(read_pairs(text),[parameters.n_frq.' real(scan).']),...
                 'RF offsets remain signed Hz coordinates rather than folded frequencies');

% Export pulse-acquire and HYSCORE arrays with native sequence parameters
parameters=struct('sweep',20e6,'offset',0,'nsteps',[3 4],'zerofill',[5 6]);
values=reshape(1:12,4,3); emr.method='HYSCORE';
text=jcamp_epr(spin_system,parameters,values,'time',emr);
result=test_true(result,'HYSCORE time grid',isequal(read_pairs(text),...
                 [repmat((0:3)'/20e6,3,1) values(:)]),...
                 'nsteps sets dimensionality and sweep sets both dwell times');
values=reshape(1:30,6,5); text=jcamp_epr(spin_system,parameters,values,'frequency',emr);
result=test_true(result,'HYSCORE frequency grid',isequal(read_pairs(text),...
                 [repmat(ft_axis(0,20e6,6).',5,1) values(:)]),...
                 'processed array lengths include zero filling without applying another transform');
parameters=struct('sweep',8e6,'offset',1e6,'npoints',5,'zerofill',5); emr.method='FID';
text=jcamp_epr(spin_system,parameters,(1:5)','time',emr);
result=test_true(result,'EPR pulse acquire',isequal(read_pairs(text),[(0:4)'/8e6 (1:5)']),...
                 'one-dimensional electron acquisition uses seconds');

% Explicit pulse-delay axes cover DEER/ESEEM and non-uniform trajectories
axes={[0;1e-8;4e-8;9e-8]}; emr.method='ELDOR';
text=jcamp_signal(spin_system,axes,scan.',{'SECOND'},{'DEER delay'},emr);
result=test_true(result,'explicit delay axis',isequal(read_pairs(text),...
                 [axes{1} real(scan).'; axes{1} imag(scan).']),...
                 'non-uniform pulse delays and phase-sensitive signals remain unchanged');

% Reject shape changes and inconsistent component dimensions before writing
parameters=struct('spins',{{'1H'}},'sweep',800,'offset',0);
result=test_throws(result,'row NMR input',@() jcamp_nmr(spin_system,parameters,1:5,{'time'},info),...
                   'one-dimensional NMR signals must be columns');
fid=struct('pos',zeros(4,3),'neg',zeros(3,4));
result=test_throws(result,'inconsistent components',@() jcamp_nmr(spin_system,parameters,fid,{'time','time'},info),...
                   'quadrature component dimensions must agree');
numeric_info=info; numeric_info.sequence=42;
result=test_throws(result,'numeric pulse description',@() jcamp_nmr(spin_system,parameters,(1:5)',{'time'},numeric_info),...
                   'pulse sequences must be described with text rather than a numeric label');
parameters.spins={'E'};
result=test_throws(result,'electron in NMR wrapper',@() jcamp_nmr(spin_system,parameters,(1:5)',{'time'},info),...
                   'electron observations belong in the EPR wrapper');
result=test_throws(result,'explicit coordinate mismatch',@() jcamp_signal(spin_system,axes,(1:5)',{'SECOND'},{'Time'},emr),...
                   'an explicit axis must match the signal length');
parameters=struct('b_axis',[0.3 0.4],'mw_freq',9e9);
result=test_throws(result,'column field-sweep input',@() jcamp_epr(spin_system,parameters,[1;2],'field',emr),...
                   'the native row scan contract is explicit rather than shape-adapted');

end

% Read uncompressed numerical pairs independently of the wrapper construction
function pairs=read_pairs(text)

% Retain only table rows, excluding headers and comments
rows=strsplit(text,char([13 10])); pairs=zeros(numel(rows),2); count=0; active=false;
for n=1:numel(rows)
    if startsWith(rows{n},'##')
        active=startsWith(rows{n},'##XYDATA=')||startsWith(rows{n},'##XYPOINTS=')||...
               startsWith(rows{n},'##DATA TABLE=');
    elseif active&&~isempty(rows{n})&&~startsWith(rows{n},'$$')
        values=sscanf(strrep(rows{n},',',' '),'%f');
        assert(numel(values)==2,'uncompressed pair row must have two values.');
        count=count+1; pairs(count,:)=values';
    end
end
pairs=pairs(1:count,:);

end

% Record a repeated refusal check without interrupting the remaining tests
function result=test_throws(result,label,operation,why)

% Call the production path and require a refusal
failed=false;
try
    operation();
catch
    failed=true;
end
result=test_true(result,label,failed,why);

end

% Consistency enforcement
function grumble()
end


