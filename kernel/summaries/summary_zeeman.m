% Prints Zeeman interaction tensor summary for a Spinach system. Syntax:
%
%                 summary_zeeman(spin_system,header)
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
% <https://spindynamics.org/wiki/index.php?title=summary_zeeman.m>

function summary_zeeman(spin_system,header)

% Check consistency
grumble(spin_system,header);

% Print the summary table
report(spin_system,header);
report(spin_system,'=======================================================================================================');
report(spin_system,'#    Spin  2S+1                      Matrix                      Tr/3         norm[rank1]  norm[rank2] ');
report(spin_system,'-------------------------------------------------------------------------------------------------------');
for n=1:spin_system.comp.nspins
    if ~isempty(spin_system.inter.zeeman.matrix{n})
        
        % Get the isotropic part
        iso=trace(spin_system.inter.zeeman.matrix{n})/3;
        
        % Get the first and second rank parts
        [~,rank1,rank2]=mat2sphten(spin_system.inter.zeeman.matrix{n});
        rank1=sphten2mat([],rank1,[]); rank2=sphten2mat([],[],rank2);
        
        % Do the printing
        report(spin_system,[pad(num2str(n),5)...
                            pad(spin_system.comp.isotopes{n},6)...
                            pad(num2str(spin_system.comp.mults(n)),7)...
                            num2str(spin_system.inter.zeeman.matrix{n}(1,:),'%+10.6e  %+10.6e  %+10.6e')]);
        report(spin_system,['                  '...
                            num2str(spin_system.inter.zeeman.matrix{n}(2,:),'%+10.6e  %+10.6e  %+10.6e')...
                            '    ' ...
                            pad(num2str(iso,'%+10.4e'),13)...
                            pad(num2str(norm(rank1,2),'%+10.4e'),13)...
                            pad(num2str(norm(rank2,2),'%+10.4e'),13)]);
        report(spin_system,['                  '...
                            num2str(spin_system.inter.zeeman.matrix{n}(3,:),'%+10.6e  %+10.6e  %+10.6e')]);
        
        % Print the break line
        if n<spin_system.comp.nspins, report(spin_system,' '); end
        
    end
end
report(spin_system,'=======================================================================================================');

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

% Jim wants me to go out with other men so that
% he will have something to write about.
%
% Nora Barnacle, the
% wife of James Joyce

