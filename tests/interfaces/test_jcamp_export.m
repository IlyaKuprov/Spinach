% Tests faithful NMR/EPR JCAMP export and refusal of ambiguous inputs.
% Syntax:
%
%                    result=test_jcamp_export()
%
% Outputs:
%
%    result - regression checks for sampling, complex channels, peak
%             assignments, linked blocks, and safe file replacement
%
% talos@spindynamics.org

function result=test_jcamp_export()

% Check consistency
grumble();

% Set up a spectrum containing values across the floating-point range
result=new_test_result('interfaces/jcamp_export','JCAMP NMR/EPR exporter',...
                       'sample coordinates, amplitudes, channels, and metadata survive export.');
block=struct('title','Test spectrum','type','NMR SPECTRUM',...
             'metadata',{{'.OBSERVE FREQUENCY',400; '.OBSERVE NUCLEUS','^1H'}},...
             'x',(5:-1:1)','y',[pi;-2;1e-200;1e200;0],...
             'xunits','HZ','yunits','ARBITRARY UNITS');
data=struct('title','Test data','origin','Spinach','owner','Public domain','blocks',{{block}});
text=jcamp_export(data);
result=test_true(result,'regular descending axis',contains(text,'##DATA CLASS=XYDATA'),...
                 'an exactly regular descending axis uses incremental data');
pairs=read_pairs(text);
result=test_true(result,'double precision',isequal(pairs,[block.x block.y]),...
                 'AFFN records recover the original IEEE double values without quantisation');
result=test_true(result,'line length',all(cellfun(@numel,strsplit(text,char([13 10])))<=80),...
                 'data and headers satisfy the JCAMP record limit');

% Integer metadata must not be rounded through a floating-point conversion
block.metadata(end+1,:)={'$COUNT',intmax('uint64')}; data.blocks={block}; text=jcamp_export(data);
result=test_true(result,'integer metadata',contains(text,'##$COUNT=18446744073709551615'),...
                 'integer metadata retains every decimal digit, including values outside the double exact-integer range');
block.metadata=block.metadata(1:2,:);

% JCAMP 5.01 permits magnitude and power labels without recomputing samples
transformed=block; transformed.y=[1;2;3;4;5]; valid=true;
for units={'MAGNITUDE','POWER'}
    transformed.yunits=units{1}; data.blocks={transformed}; text=jcamp_export(data);
    valid=valid&&contains(text,['##YUNITS=' units{1}])&&isequal(read_pairs(text),[transformed.x transformed.y]);
end
result=test_true(result,'NMR transformed units',valid,...
                 'magnitude and power ordinate units retain the already transformed caller samples');

% Irregular, non-monotonic, and singleton axes must retain explicit coordinates
block.x=[7;2;pi;-5;0]; data.blocks={block}; text=jcamp_export(data);
result=test_true(result,'irregular coordinates',contains(text,'##DATA CLASS=XYPOINTS')&&...
                 isequal(read_pairs(text),[block.x block.y]),'no interpolation, sorting, or axis fitting occurs');
block.x=pi; block.y=-2; data.blocks={block}; text=jcamp_export(data);
result=test_true(result,'singleton',contains(text,'##DATA CLASS=XYPOINTS')&&...
                 isequal(read_pairs(text),[pi -2]),'a single point has no undefined incremental spacing');

% Complex FIDs must contain two independent component pages
block.xunits='SECONDS'; block.x=(0:3)'/4; block.y=[1+2*1i;-3+4*1i;5-6*1i;-7-8*1i]; block.type='NMR FID';
block.metadata=[block.metadata; {'.DELAY','(0, 0)'; '.ACQUISITION MODE','SIMULTANEOUS'}];
data.blocks={block}; text=jcamp_export(data); pairs=read_pairs(text);
result=test_true(result,'complex channels',contains(text,'##DATA TABLE=(X++(R..R)), XYDATA')&&...
                 contains(text,'##DATA TABLE=(X++(I..I)), XYDATA')&&...
                 isequal(pairs,[block.x real(block.y); block.x imag(block.y)]),...
                 'real and imaginary amplitudes retain their signs and order');
block.x=[0;0.1;0.3;1]; data.blocks={block}; text=jcamp_export(data);
result=test_true(result,'irregular complex FID',contains(text,'##DATA TABLE=(XR..XR), PROFILE')&&...
                 isequal(read_pairs(text),[block.x real(block.y); block.x imag(block.y)]),...
                 'irregular quadrature samples use explicit pairs in NTUPLES');

% General pages represent indirect coordinates, unequal lengths, and extra channels
variables=struct('name',{'Frequency','Intensity','Time'},'symbol',{'X','Y','T'},...
                 'type',{'INDEPENDENT','DEPENDENT','INDEPENDENT'},...
                 'units',{'HZ','ARBITRARY UNITS','SECONDS'});
pages=struct('x',{[1;2;3],[1;pi;4;5]},'y',{[4;5;6],[7;8;9;10]},'xvar',{'X','X'},...
             'yvar',{'Y','Y'},'coordinates',{{'T',0},{'T',2}});
general=struct('title','General pages','type','NMR SPECTRUM',...
               'metadata',{block.metadata(1:2,:)},'variables',variables,'pages',pages);
data.blocks={general}; text=jcamp_export(data);
result=test_true(result,'ragged multidimensional pages',contains(text,'##PAGE=T=2')&&...
                 contains(text,'##NPOINTS=4')&&contains(text,'##DATA TABLE=(XY..XY), PROFILE')&&...
                 isequal(read_pairs(text),[1 4;2 5;3 6;1 7;pi 8;4 9;5 10]),...
                 'physical page coordinates and each page point count are explicit');

% Repeated coordinate sets retain separately supplied traces in page order
repeated=general; repeated.pages(2).coordinates={'T',0}; data.blocks={repeated}; text=jcamp_export(data);
result=test_true(result,'repeated page coordinates',count(text,'##PAGE=T=0')==2&&...
                 isequal(read_pairs(text),[1 4;2 5;3 6;1 7;pi 8;4 9;5 10]),...
                 'repeated coordinates do not merge, replace, or reorder independent page tables');

% Different regular grids cannot share an implicit NTUPLES abscissa
shifted=general;
shifted.pages(2).x=[10;11;12]; shifted.pages(2).y=[7;8;9];
data.blocks={shifted}; text=jcamp_export(data);
result=test_true(result,'different regular page axes',count(text,'##DATA TABLE=(XY..XY), PROFILE')==2,...
                 'page-specific endpoints are explicit when shared FIRST/LAST attributes would change them');

% Assigned peaks carry widths, multiplicity, and the peak-finding convention
peaks=struct('x',[1;2],'y',[3;4],'width',[0.1;0.2],...
             'multiplicity',{{'S';'D'}},'assignment',{{'H1';'H2'}},'method','FWHM; analytic line positions');
peak_block=struct('title','Assigned peaks','type','NMR PEAK ASSIGNMENTS',...
                  'metadata',{block.metadata(1:2,:)},'peaks',peaks,'xunits','PPM','yunits','ARBITRARY UNITS');
data.blocks={peak_block}; text=jcamp_export(data);
result=test_true(result,'assigned peaks',contains(text,'##DATA CLASS=ASSIGNMENTS')&&contains(text,'##PEAK ASSIGNMENTS=(XYMWA)')&&...
                 contains(text,'(1, 3, S, 0.10000000000000001, <H1>)')&&...
                 contains(text,'$$ FWHM; analytic line positions'),...
                 'the variable list, width convention, and assignment delimiters are retained');

% Assignment strings can contain chemical group punctuation inside angle brackets
peak_block.peaks.assignment{1}='H1(CH3), methyl'; data.blocks={peak_block}; text=jcamp_export(data);
result=test_true(result,'assignment punctuation',contains(text,'<H1(CH3), methyl>'),...
                 'parentheses and commas inside an angle-bracketed string do not become group delimiters');

% NMR assignment heights, multiplicity, and widths are independently optional
assigned=peak_block; assigned.peaks=struct('x',peaks.x,'assignment',{peaks.assignment});
data.blocks={assigned}; text=jcamp_export(data);
result=test_true(result,'position-only NMR assignments',contains(text,'##PEAK ASSIGNMENTS=(XA)'),...
                 'assigned NMR positions do not require heights or a method comment');
assigned.peaks.y=peaks.y; data.blocks={assigned}; text=jcamp_export(data);
result=test_true(result,'NMR assignments without method',contains(text,'##PEAK ASSIGNMENTS=(XYA)'),...
                 'an NMR method comment is required for widths, not ordinary assignments');
assigned.peaks.multiplicity=peaks.multiplicity; data.blocks={assigned}; text=jcamp_export(data);
result=test_true(result,'NMR multiplicity without method',contains(text,'##PEAK ASSIGNMENTS=(XYMA)'),...
                 'multiplicity does not introduce a width-convention requirement');
assigned.peaks=rmfield(assigned.peaks,'y'); assigned.peaks.width=peaks.width; assigned.peaks.method=peaks.method;
data.blocks={assigned}; text=jcamp_export(data);
result=test_true(result,'NMR width without height',contains(text,'##PEAK ASSIGNMENTS=(XMWA)'),...
                 'assigned widths and multiplicities remain optional independently of heights');

% A compound file keeps its own headers and balanced block terminators
emr=struct('title','EPR simulation','type','EMR SIMULATION',...
           'metadata',{{'.DETECTION MODE','CW'; '.METHOD','SPECTRUM';...
           '.SIMULATION SOURCE','Spinach'; '.SIMULATION PARAMETERS','Synthetic test data'}},...
           'x',[0.3;0.4],'y',[-1;1],'xunits','TESLA','yunits','ARBITRARY UNITS');
data.blocks={emr,peak_block}; text=jcamp_export(data);
result=test_true(result,'linked NMR/EPR',contains(text,'##DATA TYPE=LINK')&&contains(text,'##BLOCKS=2')&&...
                 count(text,'##TITLE=')==3&&count(text,'##END=')==3,...
                 'unrelated NMR and EPR datasets have separate complete data blocks');

% Missing observations are explicit rather than discarded or converted to zero
emr.y=[NaN;1]; data.blocks={emr}; text=jcamp_export(data);
result=test_true(result,'missing observation',contains(text,'##DATA CLASS=XYPOINTS')&&contains(text,'0.29999999999999999, ?')&&...
                 ~contains(text,'NaN'),'missing ordinates are represented by the JCAMP question mark');

% EMR peak descriptors follow the technique's single-group grammar
emr_peaks=rmfield(emr,{'x','y'}); emr_peaks.peaks=struct('x',[0.3;0.4],'y',[1;2]);
data.blocks={emr_peaks}; text=jcamp_export(data);
result=test_true(result,'EMR peak table',contains(text,'##PEAK TABLE=(XY)')&&~contains(text,'##PEAK TABLE=(XY..XY)'),...
                 'EMR peak tables use the variable list defined in the EMR protocol');
emr_peaks.peaks.assignment={'radical1';'radical2'}; emr_peaks.peaks.method='analytic';
data.blocks={emr_peaks}; text=jcamp_export(data);
result=test_true(result,'EMR assignment class',contains(text,'##DATA CLASS=PEAK ASSIGNMENTS')&&contains(text,'##PEAK ASSIGNMENTS=(XYA)'),...
                 'EMR assignments follow the explicit core-header definition and assignment grammar');

% A missing first quadrature sample must not require a finite FIRST ordinate
missing=block; missing.x=(0:3)'/4; missing.y(1)=complex(NaN,NaN);
data.blocks={missing}; text=jcamp_export(data);
result=test_true(result,'missing initial quadrature',~contains(text,')), XYDATA')&&contains(text,', PROFILE'),...
                 'leading missing samples use explicit pairs instead of incremental pages');

% Long metadata wraps into continuation records without introducing new labels
emr.metadata=[emr.metadata; {'COMMENT',repmat('sample description ',1,15)}];
data.blocks={emr}; text=jcamp_export(data);
result=test_true(result,'multiline metadata',all(cellfun(@numel,strsplit(text,char([13 10])))<=80),...
                 'long textual LDRs wrap on spaces within the record limit');

% Private metadata may use its own dotted identifier namespace
emr.metadata(end+1,:)={'$SAMPLE.DESCRIPTION','Private note'}; data.blocks={emr}; text=jcamp_export(data);
result=test_true(result,'private dotted label',contains(text,'##$SAMPLE.DESCRIPTION=Private note'),...
                 'private labels retain their user-defined namespace rather than impersonating reserved labels');

% Simulation descriptions may contain multiple lines of text
simulation=emr; simulation.metadata{3,2}={'Spinach';'Numerical simulation'};
simulation.metadata{4,2}={'Field: 0.3 T';'Linewidth: 0.001 T'};
data.blocks={simulation}; text=jcamp_export(data);
result=test_true(result,'multiline simulation descriptions',contains(text,['##.SIMULATION SOURCE=Spinach' char([13 10]) 'Numerical simulation'])&&...
                 contains(text,['##.SIMULATION PARAMETERS=Field: 0.3 T' char([13 10]) 'Linewidth: 0.001 T']),...
                 'simulation descriptions preserve explicitly supplied lines of ASCII text');

% A long value must not split the spaces inside its label
emr.metadata(end+1,:)={'SPECTROMETER/DATA SYSTEM',repmat('a',1,70)};
data.blocks={emr}; text=jcamp_export(data);
result=test_true(result,'complete label prefix',contains(text,['##SPECTROMETER/DATA SYSTEM=' char([13 10]) repmat('a',1,70)]),...
                 'an intact label precedes the long value on a continuation line');

% A requested output file must match the returned text byte for byte
file_name=[tempname() '.jdx']; data.filename=file_name;
cleanup=onCleanup(@()delete(file_name)); text=jcamp_export(data);
result=test_true(result,'file output',strcmp(fileread(file_name),text),...
                 'file output and in-memory output contain exactly the same ASCII records');
sentinel=fileread(file_name);

% A directory destination must be refused without leaving an unnamed output
folder=tempname(); mkdir(folder); folder_cleanup=onCleanup(@()rmdir(folder,'s'));
bad=data; bad.filename=folder; refused=false;
try
    jcamp_export(bad);
catch
    refused=true;
end
listing=dir(folder);
result=test_true(result,'directory destination',refused&&isempty(listing(~[listing.isdir])),...
                 'directory destinations are refused without creating an unexpected file');

% Invalid inputs must be refused before an existing destination can be modified
invalid={};
bad=data; bad.blocks{1}.y=[Inf;1]; invalid{end+1}=bad;
bad=data; bad.blocks{1}.x=[0.3 0.4]; invalid{end+1}=bad;
bad=data; bad.blocks{1}.metadata(end+1,:)={'Y_FACTOR',2}; invalid{end+1}=bad;
bad=data; bad.blocks{1}.metadata(end+1,:)={'COMMENT',sprintf('ok\n##END=')}; invalid{end+1}=bad;
bad=data; bad.blocks={general}; bad.blocks{1}.pages(2).coordinates=cell(0,2); invalid{end+1}=bad;
bad=data; bad.blocks={general}; bad.blocks{1}.variables(2).symbol='X'; invalid{end+1}=bad;
bad=data; bad.blocks={peak_block}; bad.blocks{1}.peaks.assignment{1}='<H1>'; invalid{end+1}=bad;
bad=data; bad.blocks={block}; bad.blocks{1}.metadata=cell(0,2); invalid{end+1}=bad;
bad=data; bad.blocks={block}; bad.blocks{1}.metadata{2,2}='^1Hgarbage'; invalid{end+1}=bad;
bad=data; bad.blocks={block}; bad.blocks{1}.metadata{3,2}='not-a-delay'; invalid{end+1}=bad;
bad=data; bad.blocks={block}; bad.blocks{1}.metadata{3,2}='(0, 0) trailing'; invalid{end+1}=bad;
bad=data; bad.blocks={block}; bad.blocks{1}.metadata{3,2}='(1i, 0)'; invalid{end+1}=bad;
bad=data; bad.blocks={block}; bad.blocks{1}.xunits='HZ'; invalid{end+1}=bad;
bad=data; bad.blocks={general}; bad.blocks{1}.variables(1).units='SECONDS'; invalid{end+1}=bad;
bad=data; bad.blocks={emr_peaks}; bad.blocks{1}.peaks.multiplicity={'S';'D'}; invalid{end+1}=bad;
bad=data; bad.blocks{1}.yunits='TESLA'; invalid{end+1}=bad;
bad=data; bad.blocks={general}; bad.blocks{1}.variables(2).units='TESLA'; invalid{end+1}=bad;
bad=data; bad.blocks{1}.metadata{2,2}=42; invalid{end+1}=bad;
bad=data; bad.blocks{1}.metadata{2,2}='NOT-A-METHOD'; invalid{end+1}=bad;
bad=data; bad.blocks{1}.metadata{2,2}={'SPECTRUM';'FID'}; invalid{end+1}=bad;
bad=data; bad.blocks{1}.metadata(end+1,:)={'SAMPLE.DESCRIPTION','Invalid reserved label'}; invalid{end+1}=bad;
bad=data; bad.blocks{1}.metadata{3,2}=42; invalid{end+1}=bad;
bad=data; bad.blocks{1}.metadata{4,2}=[1 2]; invalid{end+1}=bad;
bad=data; bad.blocks={block}; bad.blocks{1}.y=block.y+1i; bad.blocks{1}.yunits='MAGNITUDE'; invalid{end+1}=bad;
bad=data; bad.blocks={block}; bad.blocks{1}.y=block.y+1i; bad.blocks{1}.yunits='POWER'; invalid{end+1}=bad;
for n=1:numel(invalid)
    refused=false;
    try
        jcamp_export(invalid{n});
    catch
        refused=true;
    end
    result=test_true(result,['invalid input ' num2str(n)],refused&&strcmp(fileread(file_name),sentinel),...
                     'invalid inputs do not produce a partial file or clobber an existing destination');
end

end

% Read uncompressed pair tables independently of the exporter
function pairs=read_pairs(text)

% Only numerical rows inside a data table contribute samples
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

% Consistency enforcement
function grumble()
end


