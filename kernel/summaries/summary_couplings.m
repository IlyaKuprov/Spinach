% Prints spin-spin coupling tensor summary for a Spinach system. Syntax:
%
%                 summary_couplings(spin_system,header)
%
% Parameters:
%
%    spin_system  - Spinach spin system description object
%
%    header       - a string of text to precede the summary
%
% Outputs:
%
%    this function prints to the console or to the user-specified
%    output via report.m function
%
% ilya.kuprov@weizmann.ac.il

function summary_couplings(spin_system,header)

% Check consistency
grumble(spin_system,header);

% Print the summary table
report(spin_system,header);
report(spin_system,'======================================================================================================');
report(spin_system,'Spin A  Spin B                     Matrix                      Tr/3         norm[rank1]  norm[rank2]  ');
report(spin_system,'------------------------------------------------------------------------------------------------------');

% Detect significant couplings
[rows,cols,~]=find(cellfun(@(x)norm(x,2),spin_system.inter.coupling.matrix)>2*pi*spin_system.tols.inter_cutoff);

% Loop over significant couplings
for n=1:numel(rows)
    
    % Get the isotropic part
    iso=trace(spin_system.inter.coupling.matrix{rows(n),cols(n)})/3;
    
    % Get the first and second rank parts
    [~,rank1,rank2]=mat2sphten(spin_system.inter.coupling.matrix{rows(n),cols(n)});
    rank1=sphten2mat([],rank1,[]); rank2=sphten2mat([],[],rank2);
    
    % Do the printing in Hz
    report(spin_system,[pad(num2str(rows(n)),8)...
                        pad(num2str(cols(n)),8) ...
                        num2str(spin_system.inter.coupling.matrix{rows(n),cols(n)}(1,:)/(2*pi),'%+10.6e  %+10.6e  %+10.6e')]);
    report(spin_system,['                '...
                        num2str(spin_system.inter.coupling.matrix{rows(n),cols(n)}(2,:)/(2*pi),'%+10.6e  %+10.6e  %+10.6e')...
                        '    ' ...
                        pad(num2str(iso/(2*pi),'%+10.4e'),13)...
                        pad(num2str(norm(rank1,2)/(2*pi),'%+10.4e'),13)...
                        pad(num2str(norm(rank2,2)/(2*pi),'%+10.4e'),13)]);
    report(spin_system,['                '...
                        num2str(spin_system.inter.coupling.matrix{rows(n),cols(n)}(3,:)/(2*pi),'%+10.6e  %+10.6e  %+10.6e')]);
    
    % Print the break line
    if n<numel(rows), report(spin_system,' '); end
    
end
report(spin_system,'======================================================================================================');

end

% Consistency enforcement
function grumble(spin_system,header)
if ~isstruct(spin_system)
    error('spin_system must be a structure.');
end
if ~ischar(header)
    error('header must be a character string.');
end
end

