% Sets case-specific assumptions for various simulation contexts. This
% function determines the behaviour of the Hamiltonian generation func-
% tion and should be called before the Hamiltonian is requested. The 
% function text is self-explanatory - interaction strength parameters
% are set in each section according to the physical requirements of
% of each specific simulation context. Syntax:
%
%        spin_system=assume(spin_system,assumptions,retention)
%
% Parameters:
%
%    assumptions - 'nmr'  for high-field NMR)
%
%                  'esr'  for electron rotating frame ESR
%
%                  'deer' for DEER spectroscopy
%
%                  'deer-zz' for DEER spectroscopy with electron
%                            flip-flop terms removed
%
%                  'labframe' for full laboratory frame simulation
%                             with all Hamiltonian terms retained;
%                             bosonic modes are allowed and stay in
%                             the laboratory frame with all of their
%                             interaction terms retained
%
%                  'qnmr'     for quadrupolar NMR with numerical
%                             rotating frames: spin-1/2 particles
%                             will be in the rotating frame but
%                             spin>1/2 particles initially in the
%                             laboratory frame
%
%                  'cavity'   for cavity QED: spins and bosonic
%                             modes in a common rotating frame with
%                             the rotating wave approximation, mode
%                             energies to be built as detunings from
%                             the carrier frequency, exchange terms
%                             keeping flip-flop components only,
%                             anharmonicity, Kerr, and dispersive
%                             terms in full; longitudinal and modu-
%                             lation terms are disallowed because
%                             they average out
%
%                  'spin-phonon' for spins in their usual rotating
%                             frames with bosonic modes in the la-
%                             boratory frame: electron and nuclear
%                             terms as in the 'esr' set, electron-
%                             mode exchange terms dropped as non-
%                             secular, mode-mode and nucleus-mode
%                             exchange terms retained in full,
%                             longitudinal, dispersive, and modula-
%                             tion terms and all diagonal mode
%                             terms retained
%
%    retention  -  'zeeman' drops all spin-spin interactions 
%
%                  'couplings' drops all Zeeman interactions
%
% Outputs:
%
%    the function updates the spin_system object
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=assume.m>

function spin_system=assume(spin_system,assumptions,retention)

% Check consistency
if nargin==2
    grumble(assumptions);
elseif nargin==3
    grumble(assumptions,retention);
else
    error('incorrect number of input arguments.');
end

% Store the assumption specification
spin_system.inter.assumptions=assumptions;

% Preallocate the arrays
spin_system.inter.zeeman.strength=cell(spin_system.comp.nspins,1);
spin_system.inter.giant.strength=cell(spin_system.comp.nspins,1);
spin_system.inter.coupling.strength=cell(spin_system.comp.nspins);

% Set the approximations
switch assumptions
    
    case {'nmr'}

        % Disallow particles that are not spins
        if any(ismember({'C','V','T'},spin_system.comp.types),'all')
            error('Cavities, phonons, and transmons are not allowed under this assumption set.');
        end

        % Do the reporting
        report(spin_system,'generic high-field NMR assumption set:');
        report(spin_system,'  rotating frame approximation for electrons,');
        report(spin_system,'  rotating frame approximation for nuclei,');
        report(spin_system,'  secular terms for electron Zeeman interactions,');
        report(spin_system,'  secular terms for nuclear Zeeman interactions,');
        report(spin_system,'  secular terms for giant spin model interactions,');
        report(spin_system,'  secular coupling terms for spins belonging to the same species,');
        report(spin_system,'  weak coupling terms for spins belonging to different species,');
        report(spin_system,'  secular coupling terms for quadratic couplings.');
        
        % Process Zeeman interactions
        for n=1:spin_system.comp.nspins
            
            % All Zeeman interactions should be secular
            spin_system.inter.zeeman.strength{n}='secular';
            
        end
        
        % Process couplings
        for n=1:spin_system.comp.nspins
            
            % Giant spin terms should be secular
            spin_system.inter.giant.strength{n}='secular';
            
            % For spin-spin couplings, it depends
            for k=1:spin_system.comp.nspins
                if strcmp(spin_system.comp.isotopes{n},spin_system.comp.isotopes{k})
                    
                    % Couplings between spins of the same type should be secular
                    spin_system.inter.coupling.strength{n,k}='secular';
                    
                else
                    
                    % Couplings between spins of different types should be weak
                    spin_system.inter.coupling.strength{n,k}='zz';
                    
                end
            end
            
        end
        
    case {'esr','deer'}

        % Disallow particles that are not spins
        if any(ismember({'C','V','T'},spin_system.comp.types),'all')
            error('Cavities, phonons, and transmons are not allowed under this assumption set.');
        end
        
        % Do the reporting
        report(spin_system,'generic high-field EPR/DEER assumption set:');
        report(spin_system,'  rotating frame approximation for electrons,');
        report(spin_system,'  laboratory frame simulation for nuclei,');
        report(spin_system,'  secular terms for electron Zeeman interactions,');
        report(spin_system,'  all terms for nuclear Zeeman interactions,');
        report(spin_system,'  secular terms for inter-electron couplings,');
        report(spin_system,'  secular terms for electron zero-field splittings,');
        report(spin_system,'  secular terms for giant spin model interactions,');
        report(spin_system,'  weak and pseudosecular terms for hyperfine couplings,');
        report(spin_system,'  all terms for inter-nuclear couplings,');
        report(spin_system,'  all terms for nuclear quadrupolar couplings.');
        
        % Process Zeeman interactions
        for n=1:spin_system.comp.nspins
            if strcmp(spin_system.comp.isotopes{n}(1),'E')
                
                % Electron Zeeman interactions should be secular
                spin_system.inter.zeeman.strength{n}='secular';
                
            else
                
                % Full Zeeman tensors should be used for nuclei
                spin_system.inter.zeeman.strength{n}='full';
                
            end
        end
        
        % Process couplings
        for n=1:spin_system.comp.nspins
            
            % Giant spin terms should be secular
            spin_system.inter.giant.strength{n}='secular';
            
            % For spin-spin couplings, it depends
            for k=1:spin_system.comp.nspins
                if strcmp(spin_system.comp.isotopes{n}(1),'E')&&(~strcmp(spin_system.comp.isotopes{k}(1),'E'))
                    
                    % Couplings from electrons to nuclei should be left-secular
                    spin_system.inter.coupling.strength{n,k}='z*';
                    
                elseif strcmp(spin_system.comp.isotopes{k}(1),'E')&&(~strcmp(spin_system.comp.isotopes{n}(1),'E'))
                    
                    % Couplings from nuclei to electrons should be right-secular
                    spin_system.inter.coupling.strength{n,k}='*z';
                    
                elseif strcmp(spin_system.comp.isotopes{n}(1),'E')&&(strcmp(spin_system.comp.isotopes{k}(1),'E'))
                    
                    % Couplings between electrons should be secular
                    spin_system.inter.coupling.strength{n,k}='secular';
                    
                else
                    
                    % Couplings between nuclei should be strong
                    spin_system.inter.coupling.strength{n,k}='strong';
                    
                end
            end
        end
        
    case {'deer-zz'}

        % Disallow particles that are not spins
        if any(ismember({'C','V','T'},spin_system.comp.types),'all')
            error('Cavities, phonons, and transmons are not allowed under this assumption set.');
        end
        
        % Do the reporting
        report(spin_system,'specialised (electron ZZ only) DEER assumption set:');
        report(spin_system,'  rotating frame approximation for electrons,');
        report(spin_system,'  laboratory frame simulation for nuclei,');
        report(spin_system,'  secular terms for electron Zeeman interactions,');
        report(spin_system,'  all terms for nuclear Zeeman interactions,');
        report(spin_system,'  weak terms for inter-electron couplings,');
        report(spin_system,'  secular terms for electron zero-field splittings,');
        report(spin_system,'  secular terms for giant spin model interactions,');
        report(spin_system,'  weak and pseudosecular terms for hyperfine couplings,');
        report(spin_system,'  all terms for inter-nuclear couplings,');
        report(spin_system,'  all terms for nuclear quadrupolar couplings.');
        
        % Process Zeeman interactions
        for n=1:spin_system.comp.nspins
            if strcmp(spin_system.comp.isotopes{n}(1),'E')
                
                % Electron Zeeman interactions should be secular
                spin_system.inter.zeeman.strength{n}='secular';
                
            else
                
                % Full Zeeman tensors should be used for nuclei
                spin_system.inter.zeeman.strength{n}='full';
                
            end
        end
        
        % Process couplings
        for n=1:spin_system.comp.nspins
            
            % Giant spin terms should be secular
            spin_system.inter.giant.strength{n}='secular';
            
            % For spin-spin couplings, it depends
            for k=1:spin_system.comp.nspins
                if strcmp(spin_system.comp.isotopes{n}(1),'E')&&(~strcmp(spin_system.comp.isotopes{k}(1),'E'))
                    
                    % Couplings from electron to nucleus should be left-secular
                    spin_system.inter.coupling.strength{n,k}='z*';
                    
                elseif strcmp(spin_system.comp.isotopes{k}(1),'E')&&(~strcmp(spin_system.comp.isotopes{n}(1),'E'))
                    
                    % Couplings from nucleus to electron should be right-secular
                    spin_system.inter.coupling.strength{n,k}='*z';
                    
                elseif strcmp(spin_system.comp.isotopes{n}(1),'E')&&(strcmp(spin_system.comp.isotopes{k}(1),'E'))&&(n~=k)
                    
                    % Couplings between different electrons should be weak
                    spin_system.inter.coupling.strength{n,k}='zz';
                    
                elseif strcmp(spin_system.comp.isotopes{n}(1),'E')&&(strcmp(spin_system.comp.isotopes{k}(1),'E'))&&(n==k)
                    
                    % Quadratic couplings on each electron should be secular
                    spin_system.inter.coupling.strength{n,k}='secular';
                    
                else
                    
                    % Couplings between nuclei should be strong
                    spin_system.inter.coupling.strength{n,k}='strong';
                    
                end
            end
            
        end
        
    case {'labframe'}

        % Do the reporting
        report(spin_system,'lab frame assumption set - no assumptions:');
        report(spin_system,'  laboratory frame simulation for electrons,');
        report(spin_system,'  laboratory frame simulation for nuclei,');
        report(spin_system,'  all terms for electron Zeeman interactions,');
        report(spin_system,'  all terms for nuclear Zeeman interactions,');
        report(spin_system,'  all terms for giant spin model interactions,');
        report(spin_system,'  all terms for all couplings.');

        % Do the reporting for bosonic modes
        if any(ismember({'C','V','T'},spin_system.comp.types),'all')
            report(spin_system,'  laboratory frame simulation for bosonic modes,');
            report(spin_system,'  exchange terms keep counter-rotating components,');
            report(spin_system,'  all anharmonicity, Kerr, longitudinal, and dispersive terms retained,');
            report(spin_system,'  all Hamiltonian modulation terms retained.');
        end

        % Process Zeeman interactions
        for n=1:spin_system.comp.nspins

            % Full Zeeman tensors should be used for all spins
            spin_system.inter.zeeman.strength{n}='full';

        end

        % Process couplings
        for n=1:spin_system.comp.nspins

            % Giant spin terms should be full
            spin_system.inter.giant.strength{n}='strong';

            % The same applies to all couplings
            for k=1:spin_system.comp.nspins

                % Full coupling tensors should be used for all spins
                spin_system.inter.coupling.strength{n,k}='strong';

            end

        end

        % Process bosonic mode terms
        if isfield(spin_system.inter,'modes')
            spin_system.inter.modes.strength.frqs='full';
            spin_system.inter.modes.strength.anharms='full';
            spin_system.inter.modes.strength.exchange='strong';
            spin_system.inter.modes.strength.kerr='full';
            spin_system.inter.modes.strength.longitudinal='full';
            spin_system.inter.modes.strength.dispersive='full';
            spin_system.inter.modes.strength.coupling_mod='full';
            spin_system.inter.modes.strength.zeeman_mod='full';
        end
        
    case {'se_dnp_h+'}

        % Disallow particles that are not spins
        if any(ismember({'C','V','T'},spin_system.comp.types),'all')
            error('Cavities, phonons, and transmons are not allowed under this assumption set.');
        end
        
        % Do the reporting
        report(spin_system,'H+ part of the solid effect DNP Hamiltonian:');
        report(spin_system,'  EzNp parts of electron-nuclear couplings,');
        report(spin_system,'  T(L,+1) parts of inter-nuclear couplings,');
        report(spin_system,'  all inter-electron interactions ignored,');
        report(spin_system,'  all giant spin model interactions ignored,');
        report(spin_system,'  all quadratic interactions ignored,');
        report(spin_system,'  all Zeeman interactions ignored.');
        
        % Ignore all Zeeman interactions
        for n=1:spin_system.comp.nspins
            spin_system.inter.zeeman.strength{n}='ignore';
        end
        
        % Process couplings
        for n=1:spin_system.comp.nspins
            
            % Ignore all giant spin model terms
            spin_system.inter.giant.strength{n}='ignore';
            
            % For spin-spin couplings, it depends
            for k=1:spin_system.comp.nspins
                if n==k
                    
                    % Ignore quadratic couplings
                    spin_system.inter.coupling.strength{n,k}='ignore';

                elseif strcmp(spin_system.comp.isotopes{n}(1),'E')&&(~strcmp(spin_system.comp.isotopes{k}(1),'E'))
                    
                    % Keep EzN+ for electron-nuclear couplings
                    spin_system.inter.coupling.strength{n,k}='z+';
                    
                elseif strcmp(spin_system.comp.isotopes{k}(1),'E')&&(~strcmp(spin_system.comp.isotopes{n}(1),'E'))
                    
                    % Keep N+Ez for nuclear-electron couplings
                    spin_system.inter.coupling.strength{n,k}='+z';
                    
                elseif strcmp(spin_system.comp.isotopes{n}(1),'E')&&(strcmp(spin_system.comp.isotopes{k}(1),'E'))
                    
                    % Ignore couplings between electrons
                    spin_system.inter.coupling.strength{n,k}='ignore';
                    
                elseif (~strcmp(spin_system.comp.isotopes{k}(1),'E'))&&(~strcmp(spin_system.comp.isotopes{n}(1),'E'))
                    
                    % Keep the T(L,+1) component of inter-nuclear couplings
                    spin_system.inter.coupling.strength{n,k}='T(L,+1)';
                    
                end
            end
            
        end
        
    case {'se_dnp_h-'}

        % Disallow particles that are not spins
        if any(ismember({'C','V','T'},spin_system.comp.types),'all')
            error('Cavities, phonons, and transmons are not allowed under this assumption set.');
        end
        
        % Do the reporting
        report(spin_system,'H- part of the solid effect DNP Hamiltonian:');
        report(spin_system,'  EzNm parts of electron-nuclear couplings,');
        report(spin_system,'  T(L,-1) parts of inter-nuclear couplings,');
        report(spin_system,'  all inter-electron interactions ignored,');
        report(spin_system,'  all giant spin model interactions ignored,');
        report(spin_system,'  all quadratic interactions ignored,');
        report(spin_system,'  all Zeeman interactions ignored.');
        
        % Ignore all Zeeman interactions
        for n=1:spin_system.comp.nspins
            spin_system.inter.zeeman.strength{n}='ignore';
        end
        
        % Process couplings
        for n=1:spin_system.comp.nspins
            
            % Ignore all giant spin model terms
            spin_system.inter.giant.strength{n}='ignore';
            
            % For spin-spin couplings, it depends
            for k=1:spin_system.comp.nspins
                if n==k
                    
                    % Ignore quadratic couplings
                    spin_system.inter.coupling.strength{n,k}='ignore';

                elseif strcmp(spin_system.comp.isotopes{n}(1),'E')&&(~strcmp(spin_system.comp.isotopes{k}(1),'E'))
                    
                    % Keep EzN- for electron-nuclear couplings
                    spin_system.inter.coupling.strength{n,k}='z-';
                    
                elseif strcmp(spin_system.comp.isotopes{k}(1),'E')&&(~strcmp(spin_system.comp.isotopes{n}(1),'E'))
                    
                    % Keep N-Ez for nuclear-electron couplings
                    spin_system.inter.coupling.strength{n,k}='-z';
                    
                elseif strcmp(spin_system.comp.isotopes{n}(1),'E')&&(strcmp(spin_system.comp.isotopes{k}(1),'E'))
                    
                    % Ignore couplings between electrons
                    spin_system.inter.coupling.strength{n,k}='ignore';
                    
                elseif (~strcmp(spin_system.comp.isotopes{k}(1),'E'))&&(~strcmp(spin_system.comp.isotopes{n}(1),'E'))
                    
                    % Keep the T(L,-1) component of inter-nuclear couplings
                    spin_system.inter.coupling.strength{n,k}='T(L,-1)';
                    
                end
            end
            
        end
        
    case {'se_dnp_h0'}

        % Disallow particles that are not spins
        if any(ismember({'C','V','T'},spin_system.comp.types),'all')
            error('Cavities, phonons, and transmons are not allowed under this assumption set.');
        end
        
        % Do the reporting
        report(spin_system,'H0 part of the solid effect DNP Hamiltonian:');
        report(spin_system,'  EzNz parts of electron-nuclear couplings,');
        report(spin_system,'  secular parts of inter-nuclear couplings,');
        report(spin_system,'  all inter-electron interactions ignored,');
        report(spin_system,'  all giant spin model interactions ignored,');
        report(spin_system,'  all quadratic interactions ignored,');
        report(spin_system,'  all Zeeman interactions ignored.');
        
        % Ignore all Zeeman interactions
        for n=1:spin_system.comp.nspins
            spin_system.inter.zeeman.strength{n}='ignore';
        end
        
        % Process couplings
        for n=1:spin_system.comp.nspins
            
            % Ignore all giant spin model terms
            spin_system.inter.giant.strength{n}='ignore';
            
            % For spin-spin couplings, it depends
            for k=1:spin_system.comp.nspins
                if n==k
                    
                    % Ignore quadratic couplings
                    spin_system.inter.coupling.strength{n,k}='ignore';

                elseif strcmp(spin_system.comp.isotopes{n}(1),'E')&&(~strcmp(spin_system.comp.isotopes{k}(1),'E'))
                    
                    % Keep EzNz for electron-nuclear couplings
                    spin_system.inter.coupling.strength{n,k}='zz';
                    
                elseif strcmp(spin_system.comp.isotopes{k}(1),'E')&&(~strcmp(spin_system.comp.isotopes{n}(1),'E'))
                    
                    % Keep NzEz for nuclear-electron couplings
                    spin_system.inter.coupling.strength{n,k}='zz';
                    
                elseif strcmp(spin_system.comp.isotopes{n}(1),'E')&&(strcmp(spin_system.comp.isotopes{k}(1),'E'))
                    
                    % Ignore couplings between electrons
                    spin_system.inter.coupling.strength{n,k}='ignore';
                    
                elseif (~strcmp(spin_system.comp.isotopes{k}(1),'E'))&&(~strcmp(spin_system.comp.isotopes{n}(1),'E'))
                    
                    % Keep the secular component of inter-nuclear couplings
                    spin_system.inter.coupling.strength{n,k}='secular';
                    
                end
            end
            
        end
        
    case 'qnmr'

        % Disallow particles that are not spins
        if any(ismember({'C','V','T'},spin_system.comp.types),'all')
            error('Cavities, phonons, and transmons are not allowed under this assumption set.');
        end
        
        % Do the reporting
        report(spin_system,'generic quadrupolar NMR assumption set:');
        report(spin_system,'  rotating frame approximation for S=1/2 nuclei,');
        report(spin_system,'  laboratory frame simulation for S>1/2 nuclei,');
        report(spin_system,'  electrons disallowed by the approximation,');
        report(spin_system,'  secular Zeeman terms for S=1/2 nuclei,');
        report(spin_system,'  all Zeeman terms for S>1/2 nuclei,');
        report(spin_system,'  secular terms for couplings between S=1/2 nuclei,');
        report(spin_system,'  weak and pseudosecular terms for couplings between S=1/2 and S>1/2 nuclei,');
        report(spin_system,'  all terms for couplings between S>1/2 nuclei,');
        report(spin_system,'  all terms for nuclear quadrupolar couplings.');
        
        % Process Zeeman interactions
        for n=1:spin_system.comp.nspins
            
            % Check for electrons
            if strcmp(spin_system.comp.isotopes{n}(1),'E')
                error('electrons are not permitted in overtone NMR simuilations.');
            end
            
            % Assign strength parameters
            if spin_system.comp.mults(n)==2
                
                % S=1/2 Zeeman interactions should be secular
                spin_system.inter.zeeman.strength{n}='secular';
                
            else
                
                % All other Zeeman interactions should be strong
                spin_system.inter.zeeman.strength{n}='full';
                
            end
            
        end
        
        % Process couplings
        for n=1:spin_system.comp.nspins
            
            % Ignore all giant spin model terms
            spin_system.inter.giant.strength{n}='ignore';
            
            % For spin-spin couplings, it depends
            for k=1:spin_system.comp.nspins
                
                % Assign strength parameters
                if (spin_system.comp.mults(n)==2)&&(spin_system.comp.mults(k)>2)
                    
                    % Couplings from S=1/2 to S>1/2 nuclei should be left-secular
                    spin_system.inter.coupling.strength{n,k}='z*';
                    
                elseif (spin_system.comp.mults(n)>2)&&(spin_system.comp.mults(k)==2)
                    
                    % Couplings from S>1/2 to S=1/2 nuclei should be right-secular
                    spin_system.inter.coupling.strength{n,k}='*z';
                    
                elseif (spin_system.comp.mults(n)==2)&&(spin_system.comp.mults(k)==2)
                    
                    % Couplings between S=1/2 nuclei should be secular
                    spin_system.inter.coupling.strength{n,k}='secular';
                    
                else
                    
                    % Couplings between S>1/2 nuclei should be strong
                    spin_system.inter.coupling.strength{n,k}='strong';
                    
                end
                
            end
            
        end

    case {'cavity'}

        % Do the reporting
        report(spin_system,'cavity QED assumption set - common rotating frame with RWA:');
        report(spin_system,'  rotating frame approximation for all spins,');
        report(spin_system,'  secular terms for all Zeeman interactions,');
        report(spin_system,'  secular terms for giant spin model interactions,');
        report(spin_system,'  secular terms for all spin-spin couplings,');
        report(spin_system,'  mode energies to be built as carrier detunings,');
        report(spin_system,'  exchange terms keep flip-flop components only,');
        report(spin_system,'  anharmonicity, Kerr, and dispersive terms retained in full.');

        % Process Zeeman interactions
        for n=1:spin_system.comp.nspins

            % All Zeeman interactions should be secular
            spin_system.inter.zeeman.strength{n}='secular';

        end

        % Process couplings
        for n=1:spin_system.comp.nspins

            % Giant spin terms should be secular
            spin_system.inter.giant.strength{n}='secular';

            % All spin-spin couplings should be secular
            for k=1:spin_system.comp.nspins
                spin_system.inter.coupling.strength{n,k}='secular';
            end

        end

        % Process bosonic mode terms
        if isfield(spin_system.inter,'modes')

            % Disallow longitudinal terms that average out under the RWA
            if ~all(cellfun(@isempty,spin_system.inter.modes.longitudinal(:)))
                error('longitudinal couplings average to zero in this frame, use spin-phonon or labframe.');
            end

            % Disallow modulation terms that average out under the RWA
            if (~all(cellfun(@isempty,spin_system.inter.modes.coupling_mod(:))))||...
               (~all(cellfun(@isempty,spin_system.inter.modes.zeeman_mod(:))))
                error('Hamiltonian modulation terms average to zero in this frame, use spin-phonon or labframe.');
            end

            % Set mode term strengths
            spin_system.inter.modes.strength.frqs='offset';
            spin_system.inter.modes.strength.anharms='full';
            spin_system.inter.modes.strength.exchange='rwa';
            spin_system.inter.modes.strength.kerr='full';
            spin_system.inter.modes.strength.longitudinal='ignore';
            spin_system.inter.modes.strength.dispersive='full';
            spin_system.inter.modes.strength.coupling_mod='ignore';
            spin_system.inter.modes.strength.zeeman_mod='ignore';

        end

    case {'spin-phonon'}

        % Do the reporting
        report(spin_system,'spin-phonon assumption set - lab frame bosonic modes:');
        report(spin_system,'  rotating frame approximation for electrons,');
        report(spin_system,'  laboratory frame simulation for nuclei,');
        report(spin_system,'  laboratory frame simulation for bosonic modes,');
        report(spin_system,'  secular terms for electron Zeeman interactions,');
        report(spin_system,'  all terms for nuclear Zeeman interactions,');
        report(spin_system,'  secular terms for inter-electron couplings,');
        report(spin_system,'  secular terms for giant spin model interactions,');
        report(spin_system,'  weak and pseudosecular terms for hyperfine couplings,');
        report(spin_system,'  all terms for inter-nuclear couplings,');
        report(spin_system,'  electron-mode exchange terms dropped as non-secular in this frame,');
        report(spin_system,'  mode-mode and nucleus-mode exchange terms retained in full,');
        report(spin_system,'  all longitudinal, dispersive, and modulation terms retained.');

        % Process Zeeman interactions
        for n=1:spin_system.comp.nspins
            if strcmp(spin_system.comp.isotopes{n}(1),'E')

                % Electron Zeeman interactions should be secular
                spin_system.inter.zeeman.strength{n}='secular';

            else

                % Full Zeeman tensors should be used for nuclei
                spin_system.inter.zeeman.strength{n}='full';

            end
        end

        % Process couplings
        for n=1:spin_system.comp.nspins

            % Giant spin terms should be secular
            spin_system.inter.giant.strength{n}='secular';

            % For spin-spin couplings, it depends
            for k=1:spin_system.comp.nspins
                if strcmp(spin_system.comp.isotopes{n}(1),'E')&&(~strcmp(spin_system.comp.isotopes{k}(1),'E'))

                    % Couplings from electrons to nuclei should be left-secular
                    spin_system.inter.coupling.strength{n,k}='z*';

                elseif strcmp(spin_system.comp.isotopes{k}(1),'E')&&(~strcmp(spin_system.comp.isotopes{n}(1),'E'))

                    % Couplings from nuclei to electrons should be right-secular
                    spin_system.inter.coupling.strength{n,k}='*z';

                elseif strcmp(spin_system.comp.isotopes{n}(1),'E')&&(strcmp(spin_system.comp.isotopes{k}(1),'E'))

                    % Couplings between electrons should be secular
                    spin_system.inter.coupling.strength{n,k}='secular';

                else

                    % Couplings between nuclei should be strong
                    spin_system.inter.coupling.strength{n,k}='strong';

                end
            end
        end

        % Process bosonic mode terms
        if isfield(spin_system.inter,'modes')
            spin_system.inter.modes.strength.frqs='full';
            spin_system.inter.modes.strength.anharms='full';
            spin_system.inter.modes.strength.exchange='nonelec';
            spin_system.inter.modes.strength.kerr='full';
            spin_system.inter.modes.strength.longitudinal='full';
            spin_system.inter.modes.strength.dispersive='full';
            spin_system.inter.modes.strength.coupling_mod='full';
            spin_system.inter.modes.strength.zeeman_mod='full';
        end

    otherwise
        
        % Complain and bomb out
        error('unrecognized assumption set.')
        
end

% Refuse retention options when bosonic modes are present
if exist('retention','var')&&ismember(retention,{'zeeman','couplings'})&&...
   isfield(spin_system.inter,'modes')
    error('retention options are undefined for spin-boson systems.');
end

% Choose the terms to retain
if exist('retention','var')&&strcmp(retention,'couplings')
    for n=1:spin_system.comp.nspins
        spin_system.inter.zeeman.strength{n}='ignore';
    end
    report(spin_system,'WARNING - all Zeeman interactions will be ignored.');
elseif exist('retention','var')&&strcmp(retention,'zeeman')
    for n=1:spin_system.comp.nspins
        spin_system.inter.giant.strength{n}='ignore';
        for k=1:spin_system.comp.nspins
            spin_system.inter.coupling.strength{n,k}='ignore';
        end
    end
    report(spin_system,'WARNING - only Zeeman interactions will be included.');
end

end

% Consistency enforcement
function grumble(assumptions,retention)
if ~ischar(assumptions)
    error('assumptions argument muct be a character string.');
end
if exist('retention','var')&&(~ischar(retention))
    error('retention argument must be a character string.');
end
end

% Humanity has advanced, when it has advanced, not because it has
% been sober, responsible, and cautious, but because it has been
% playful, rebellious, and immature.
%
% Tom Robbins

