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
%
% <https://spindynamics.org/wiki/index.php?title=summary_couplings.m>

function summary_couplings(spin_system,header)

% Check consistency
grumble(spin_system,header);
if strcmp(spin_system.sys.output,'hush'), return; end

% Print the summary table
report(spin_system,header);
report(spin_system,'======================================================================================================');
report(spin_system,'Spin A  Spin B                     Matrix                      Tr/3         norm[rank1]  norm[rank2]  ');
report(spin_system,'------------------------------------------------------------------------------------------------------');

% Detect significant couplings
pairs=spin_system.inter.coupling.matrix;
pairs=pairs(cellfun(@(x)norm(x,2),pairs(:,3))>2*pi*spin_system.tols.inter_cutoff,:);

% Loop over significant couplings
for n=1:size(pairs,1)
    
    % Get the isotropic part
    iso=trace(pairs{n,3})/3;
    
    % Get the first and second rank parts
    [~,rank1,rank2]=mat2sphten(pairs{n,3});
    rank1=sphten2mat([],rank1,[]); rank2=sphten2mat([],[],rank2);
    
    % Do the printing in Hz
    report(spin_system,[pad(num2str(pairs{n,1}),8)...
                        pad(num2str(pairs{n,2}),8) ...
                        num2str(pairs{n,3}(1,:)/(2*pi),'%+10.6e  %+10.6e  %+10.6e')]);
    report(spin_system,['                '...
                        num2str(pairs{n,3}(2,:)/(2*pi),'%+10.6e  %+10.6e  %+10.6e')...
                        '    ' ...
                        pad(num2str(iso/(2*pi),'%+10.4e'),13)...
                        pad(num2str(norm(rank1,2)/(2*pi),'%+10.4e'),13)...
                        pad(num2str(norm(rank2,2)/(2*pi),'%+10.4e'),13)]);
    report(spin_system,['                '...
                        num2str(pairs{n,3}(3,:)/(2*pi),'%+10.6e  %+10.6e  %+10.6e')]);
    
    % Print the break line
    if n<size(pairs,1), report(spin_system,' '); end
    
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

% A gentleman a few rows in front of us took grave
% exception to the behaviour of an opposing player
% and identified him, very loudly, as the author
% of The Critique of Pure Reason - repeatedly and
% with venom.
%
% Toby Young


