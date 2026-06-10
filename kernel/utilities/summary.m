% Prints various summaries on behalf of the spin system setup modules
% of Spinach kernel. Only the parameters supplied by the user are re-
% turned (if supplied); nothing new is calculated. Syntax:
%
%                  summary(spin_system,topic,header)
%
% Parameters:
%
%     spin_system   - the spin system object for which the 
%                     summary is needed
%
%     topic         - topic of the summary: 'zeeman',
%                     'pbc', 'couplings', 'rlx_rates_t1_t2', 
%                     'rlx_rates_lindblad', 'rlx_rates_nott',
%                     'rlx_rates_weiz', 'symmetry', 'basis',
%                     'basis_settings', 'coordinates'
%
%     header        - a string of text to precede the summary
%
% Outputs:
%
%     this function prints to the console or to the user specified
%     output via report.m function
%
% Note: all output produced by this function may be silenced by set-
%       ting sys.output='hush' in the Spinach input stream or by chan-
%       ging spin_system.sys.output='hush' at any point during the
%       calculation.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=summary.m>

function summary(spin_system,topic,header)

% The default header is empty
if ~exist('header','var'), header=''; end

% Check consistency
grumble(topic,header);

% If console output is disabled, simply exit
if strcmp(spin_system.sys.output,'hush'), return; end

% Do the necessary
switch topic
    
    case 'zeeman'
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
        
    case 'coordinates'
        report(spin_system,header);
        report(spin_system,'======================================');
        report(spin_system,'N    Spin     X         Y         Z   ');
        report(spin_system,'--------------------------------------');
        for n=1:spin_system.comp.nspins
            report(spin_system,[strjust([num2str(n) blanks(3-length(num2str(n)))],'left') ' '...
                                strjust([spin_system.comp.isotopes{n} blanks(5-length(spin_system.comp.isotopes{n}))],'center') '  '...
                                num2str(spin_system.inter.coordinates{n},'%+5.3f    ') '  ' spin_system.comp.labels{n}]);
        end
        report(spin_system,'======================================');
        
    case 'pbc'
        report(spin_system,header);
        report(spin_system,'===============================');
        report(spin_system,'     X         Y         Z     ');
        report(spin_system,'-------------------------------');
        for n=1:numel(spin_system.inter.pbc)
            report(spin_system,['   ' pad(num2str(spin_system.inter.pbc{n}(1),'%+5.3f   '),10)...
                                      pad(num2str(spin_system.inter.pbc{n}(2),'%+5.3f   '),10)...
                                      pad(num2str(spin_system.inter.pbc{n}(3),'%+5.3f   '),10)]);
        end
        report(spin_system,'===============================');
        
    case 'couplings'
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
         
    case 'chemistry'

        % Report multiple chemical subsystems
        if numel(spin_system.chem.parts)>1

            % Report spin system partitioning
            for n=1:numel(spin_system.chem.parts)
                report(spin_system,['chemical subsystem ' num2str(n) ' contains spins: ' num2str(spin_system.chem.parts{n})]);
            end

            % Report first-order reaction rates
            if isfield(spin_system.chem,'rates')
                report(spin_system,'inter-subsystem reaction rates:');
                report(spin_system,'===============================');
                report(spin_system,' N(from)   N(to)    Rate(Hz)   ');
                report(spin_system,'-------------------------------');
                [rows,cols,vals]=find(spin_system.chem.rates);
                for n=1:length(vals)
                    report(spin_system,[' ' strjust([num2str(rows(n)) blanks(3-length(num2str(rows(n))))],'left') '       '...
                                            strjust([num2str(cols(n)) blanks(3-length(num2str(cols(n))))],'left') '      '...
                                                     num2str(vals(n),'%+0.3e')]);
                end
                report(spin_system,'===============================');
            end

        end

        % Report flux rates if specified
        if isfield(spin_system.chem,'flux_rate')
            [rows,cols,vals]=find(spin_system.chem.flux_rate);
            if numel(vals)>0
                report(spin_system,'point-to-point flux rates:');
                report(spin_system,'===============================');
                report(spin_system,' N(from)   N(to)    Rate(Hz)   ');
                report(spin_system,'-------------------------------');
                for n=1:length(vals)
                    report(spin_system,[' ' strjust([num2str(rows(n)) blanks(3-length(num2str(rows(n))))],'left') '       '...
                                            strjust([num2str(cols(n)) blanks(3-length(num2str(cols(n))))],'left') '      '...
                                                     num2str(vals(n),'%+0.3e')]);
                end
            end
        end

    case 'rlx_rates_t1_t2'
        report(spin_system,header);
        report(spin_system,'========================================');
        report(spin_system,'N    Spin        R1             R2      ');
        report(spin_system,'----------------------------------------');
        for n=1:spin_system.comp.nspins
            report(spin_system,[strjust([num2str(n) blanks(3-length(num2str(n)))],'left') ' '...
                                strjust([spin_system.comp.isotopes{n} blanks(5-length(spin_system.comp.isotopes{n}))],'center') '  '...
                                num2str(spin_system.rlx.r1_rates(n),'%+0.5e   ') '  '...
                                num2str(spin_system.rlx.r2_rates(n),'%+0.5e   ') '  '...
                                spin_system.comp.labels{n}]);
        end
        report(spin_system,'========================================');
        
    case 'rlx_rates_lindblad'
        report(spin_system,header);
        report(spin_system,'========================================');
        report(spin_system,'N    Spin        R1             R2      ');
        report(spin_system,'----------------------------------------');
        for n=1:spin_system.comp.nspins
            report(spin_system,[strjust([num2str(n) blanks(3-length(num2str(n)))],'left') ' '...
                                strjust([spin_system.comp.isotopes{n} blanks(5-length(spin_system.comp.isotopes{n}))],'center') '  '...
                                num2str(spin_system.rlx.lind_r1_rates(n),'%+0.5e   ') '  '...
                                num2str(spin_system.rlx.lind_r2_rates(n),'%+0.5e   ') '  '...
                                spin_system.comp.labels{n}]);
        end
        report(spin_system,'========================================');
        
    case 'rlx_rates_nott'
        report(spin_system,' ');
        report(spin_system,'=== Nottingham DNP relaxation theory ===');
        report(spin_system,['Electron R1: ' num2str(spin_system.rlx.nott_r1e) ' Hz']);
        report(spin_system,['Electron R2: ' num2str(spin_system.rlx.nott_r2e) ' Hz']);
        report(spin_system,['Nuclear R1:  ' num2str(spin_system.rlx.nott_r1n) ' Hz']);
        report(spin_system,['Nuclear R2:  ' num2str(spin_system.rlx.nott_r2n) ' Hz']);
        report(spin_system,'========================================');
        
    case 'rlx_rates_weiz'
        report(spin_system,' ');
        report(spin_system,'==== Weizmann DNP relaxation theory ====');
        report(spin_system,['Electron R1: ' num2str(spin_system.rlx.weiz_r1e) ' Hz']);
        report(spin_system,['Electron R2: ' num2str(spin_system.rlx.weiz_r2e) ' Hz']);
        report(spin_system,['Nuclear R1:  ' num2str(spin_system.rlx.weiz_r1n) ' Hz']);
        report(spin_system,['Nuclear R2:  ' num2str(spin_system.rlx.weiz_r2n) ' Hz']);
        [rows,cols,vals]=find(spin_system.rlx.weiz_r1d);
        for n=1:numel(vals)
            report(spin_system,['Inter-nuclear dipolar R1(' num2str(rows(n)) ',' num2str(cols(n)) '): ' num2str(vals(n)) ' Hz']);
        end
        [rows,cols,vals]=find(spin_system.rlx.weiz_r2d);
        for n=1:numel(vals)
            report(spin_system,['Inter-nuclear dipolar R2(' num2str(rows(n)) ',' num2str(cols(n)) '): ' num2str(vals(n)) ' Hz']);
        end
        report(spin_system,'========================================');
        
    case 'symmetry'
        report(spin_system,header);
        report(spin_system,'=====================');
        report(spin_system,' Group    Spins      ');
        report(spin_system,'---------------------');
        for n=1:length(spin_system.bas.sym_spins)
            report(spin_system,['  ' spin_system.bas.sym_group{n} '     ' ...
                                     num2str(spin_system.bas.sym_spins{n})]);
        end
        report(spin_system,'====================='); report(spin_system,' ');
        
    case 'basis_summary' % Detailed summary of spherical tensor basis sets

        % Loop over substances
        for n=1:spin_system.bas.nsubst
            
            % Header for the current substance
            report(spin_system, '=============================================================');
            report(spin_system,[' Basis set, substance ' int2str(n) ': (L,M) in IST products']);
            report(spin_system, '-------------------------------------------------------------');

            % Decide if the basis should be printed
            if spin_system.bas.nstates(n)>spin_system.tols.basis_hush

                % Report that the basis set table is too long
                report(spin_system,[' over ' num2str(spin_system.tols.basis_hush) ...
                                    ' states in the basis - printing suppressed.']);

            else
     
                % Get spin list and count for current substance
                spins_in_subst=spin_system.chem.parts{n}(:);
                nspins_in_subst=numel(spins_in_subst);

                % Get single-particle ranks and projections
                [sp_ranks,sp_projs]=lin2lm(spin_system.bas.basis{n});

                % Single-particle rank range
                max_sp_rank=max(sp_ranks,[],'all');
                min_sp_rank=min(sp_ranks,[],'all');
        
                % Single-particle proj range
                max_sp_proj=max(sp_projs,[],'all');
                min_sp_proj=min(sp_projs,[],'all');
        
                % Decide basis specification table column width
                MxR=int2str(max_sp_rank); MnR=int2str(min_sp_rank);
                MxP=int2str(max_sp_proj); MnP=int2str(min_sp_proj);
                colw=1+1+max([numel(MxR) numel(MnR)])+...
                       1+max([numel(MxP) numel(MnP)])+1+1;

                % Construct and print the header line with spin numbers
                header_line=['N' repmat(blanks(colw),1,nspins_in_subst-1)];
                for k=1:nspins_in_subst
                    spin_number=int2str(spins_in_subst(k));
                    l_pos=colw*k-numel(spin_number)+1; r_pos=colw*k;
                    header_line(l_pos:r_pos)=spin_number;
                end
                report(spin_system,header_line);
                
                % Over states in the current substance
                for k=1:spin_system.bas.nstates(n)

                    % Preallocate the state specification line
                    current_line=repmat(blanks(colw),1,nspins_in_subst);
                 
                    % Over spins in the current substance
                    for m=1:nspins_in_subst

                        % Assemble the current state specification 
                        current_spec=['(' int2str(sp_ranks(k,m)) ...
                                      ',' int2str(sp_projs(k,m)) ')'];

                        % Write the state specification string
                        l_pos=colw*m-numel(current_spec)+1; r_pos=colw*m;
                        current_line(l_pos:r_pos)=current_spec;
                  
                    end

                    % Print the state specification string
                    report(spin_system,[blanks(floor((colw-4)/2)) current_line]);

                end

            end

            % Terminate the current substance basis set table
            report(spin_system,'-------------------------------------------------------------');
            report(spin_system,' ');

        end
        
    case 'basis_settings'
        
        % Report formalism and layout
        switch spin_system.bas.formalism

            case 'zeeman-wavef'

                % Schrodinger equation for the wavefunction
                report(spin_system,'Zeeman basis set using wavefunction formalism:');
                report(spin_system,'operators are matrices, states are vectors.');

            case 'zeeman-hilb'

                % Density operator formalism with LvN equation
                report(spin_system,'Zeeman basis set using Hilbert space formalism:');
                report(spin_system,'operators are matrices, states are matrices.');

            case 'zeeman-liouv'

                % Adjoint representation of the density operator formalism
                report(spin_system,'Zeeman basis set using Liouville space formalism:');
                report(spin_system,'superoperators are matrices, states are vectors.');

            case 'sphten-liouv'

                % Adjoint representation of the density oeprator formalism, IST basis set
                report(spin_system,'Spherical tensor basis set using Liouville space formalism:');
                report(spin_system,'superoperators are matrices, states are vectors.');

            otherwise

                % Complain and bomb out
                error('unrecognised formalism - see the basis section of the manual.');

        end
        
        % Report basis set approximations
        if strcmp(spin_system.bas.formalism,'sphten-liouv')
            switch spin_system.bas.approximation
                case 'IK-0'
                    report(spin_system,['IK-0 approximation | All product states between up to ' int2str(spin_system.bas.inter_level) ' spins']);
                    report(spin_system, 'IK-0 approximation | within the same chemical substance, irrespec-');
                    report(spin_system, 'IK-0 approximation | tive of proximity or interaction amplitude.');
                case 'IK-1'
                    report(spin_system,['IK-1 approximation | Product states between up to ' int2str(spin_system.bas.inter_level) ' interacting ']);
                    report(spin_system, 'IK-1 approximation | spins; interaction graph to be obtained using');
                    if strcmp(spin_system.bas.connectivity,'scalar_couplings')
                        report(spin_system,'IK-1 approximation | isotropic parts of interaction tensors.');
                    elseif strcmp(spin_system.bas.connectivity,'full_tensors')
                        report(spin_system,'IK-1 approximation | full 3x3 interaction tensors.');
                    else 
                        error('unrecognised connectivity type - see basis set preparation manual.');
                    end
                    report(spin_system, 'IK-1 approximation |');
                    report(spin_system,['IK-1 approximation | Interaction tensor drop tolerance (2-norm): ' num2str(spin_system.tols.inter_cutoff) ' Hz']);
                    report(spin_system, 'IK-1 approximation |');
                    report(spin_system,['IK-1 approximation | Product states between up to ' int2str(spin_system.bas.prox_level) ' spins']);
                    report(spin_system,['IK-1 approximation | within ' num2str(spin_system.tols.prox_cutoff) ' Angstrom of each other.']);
                case 'IK-2'
                    report(spin_system, 'IK-2 approximation | Product states involving interaction partners of');
                    report(spin_system, 'IK-2 approximation | each spin; interaction graph to be obtained using');
                    if strcmp(spin_system.bas.connectivity,'scalar_couplings')
                        report(spin_system,'IK-2 approximation | isotropic parts of interaction tensors.');
                    elseif strcmp(spin_system.bas.connectivity,'full_tensors')
                        report(spin_system,'IK-2 approximation | full 3x3 interaction tensors.');
                    else 
                        error('unrecognised connectivity type - see basis set preparation manual.');
                    end
                    report(spin_system, 'IK-2 approximation |');
                    report(spin_system,['IK-2 approximation | Interaction tensor drop tolerance (2-norm): ' num2str(spin_system.tols.inter_cutoff) ' Hz']);
                    report(spin_system, 'IK-2 approximation |');
                    report(spin_system,['IK-2 approximation | Product states between up to ' int2str(spin_system.bas.prox_level) ' spins within']);
                    report(spin_system,['IK-2 approximation | ' num2str(spin_system.tols.prox_cutoff) ' Angstrom of each other.']);
                case 'IK-DNP'
                    report(spin_system,['IK-DNP approximation | max inter-electron correlation level:   ' int2str(spin_system.bas.inter_level(1))]);
                    report(spin_system,['IK-DNP approximation | max electron-nuclear correlation level: ' int2str(spin_system.bas.inter_level(2))]);
                    report(spin_system,['IK-DNP approximation | max inter-nuclear correlation level:    ' int2str(spin_system.bas.inter_level(3))]);
                    report(spin_system, 'IK-DNP approximation | with nearest neighbours on the interaction graph.');
                    report(spin_system, 'IK-DNP approximation |');
                    report(spin_system,['IK-DNP approximation | Interaction tensor drop tolerance (2-norm): ' num2str(spin_system.tols.inter_cutoff) ' Hz']);
                case 'none'
                    report(spin_system, 'starting with the complete basis set...');
                otherwise
                    error('unrecognized approximation level - see basis set preparation manual.');
            end
        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown topic.');
        
end

end

% Consistency enforcement
function grumble(topic,header)
if ~ischar(topic)
    error('topic must be a character string.');
end
if ~ischar(header)
    error('header must be a character string.');
end
end

% 1: Blessed are the strong, for they shall possess the earth -- cursed are
% the weak, for they shall inherit the yoke.
%
% 2: Blessed are the powerful, for they shall be reverenced among men --
% cursed are the feeble, for they shall be blotted out.
%
% 3: Blessed are the bold, for they shall be masters of the world -- cursed
% are the humble, for they shall be trodden under hoofs.
%
% 4: Blessed are the victorious, for victory is the basis of right --
% cursed are the vanquished, for they shall be vassals forever.
%
% 5: Blessed are the iron-handed, for the unfit shall flee before them --
% cursed are the poor in spirit, for they shall be spat upon.
%
% 6: Blessed are the death-defiant, for their days shall be long in the
% lands -- cursed are the gazers toward a richer life beyond the grave, for
% they shall perish amidst plenty.
%
% 7: Blessed are the destroyers of false hope, for they are true Messiahs --
% cursed are the God-adorers, for they shall be shorn sheep.
%
% 8: Blessed are the valiant, for they shall obtain great treasure --
% cursed are the believers in good and evil, for they are frightened by
% shadows.
%
% 9: Blessed are those who believe in what is best for them, for never
% shall their minds be terrorized -- cursed are the "lambs of God", for
% they shall be bled whiter than snow.
%
% 10: Blessed is the man who has a sprinkling of enemies, for they shall
% make him a hero -- cursed is he who doeth good unto others who sneer upon
% him in return, for he shall be despised.
%
% 11: Blessed are the mighty-minded, for they shall ride the whirlwinds --
% cursed are they who teach lies for truth and truth for lies, for they are
% an abomination.
%
% 12: Thrice cursed are the weak whose insecurity makes them vile, for they
% shall serve and suffer.
%
% 13: The angel of self-deceit is camped in the souls of the "righteous" --
% the eternal flame of power through joy dwelleth within the flesh of a
% Satanist.
%
% Anton Szandor LaVey, "Satanic Bible"

