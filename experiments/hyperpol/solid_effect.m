% Solid effect DNP experiment, computed using the large-scale formalism
% described in (http://dx.doi.org/10.1039/C2CP23233B). The system is re-
% stricted to one electron and one nucleus type, but the number of nuc-
% lei may be very large. Syntax:
%
%               answer=solid_effect(spin_system,parameters)
%  
% Parameters:
%
%           parameters.mw_pwr   - microwave power in rad/s
%
%           parameters.theory   - level of theory. Set to 'exact' for 
%                                 the electron rotating frame calcula-
%                                 tion or to any of the following six 
%                                 options for the average Hamiltonian
%                                 theory calculation on top of the el-
%                                 ectron + nuclear rotating frame:
%                                 'ah_first_order', 'ah_second_order',
%                                 'ah_third_order', 'kb_first_order',
%                                 'kb_second_order', 'kb_third_order'.
%                                 See average.m function for the mea-
%                                 ning of these options.
%
%      parameters.nuclear_frq   - nuclear Zeeman frequency in rad/s
%
%        parameters.calc_type   - set to 'time_dependence' to get the 
%                                 time dependence of the longitudinal
%                                 magnetization and to 'steady_state'
%                                 to get the asymptotic longitudinal
%                                 magnetization.
%
%        parameters.time_step   - if 'time_dependence' is set in the
%                                 calc_type parameter, sets the time
%                                 step, seconds.
%
%          parameters.n_steps   - if 'time_dependence' is set in the
%                                 calc_type parameter, sets the num-
%                                 ber of time steps.
%
% Outputs:
%
%      answer  -  with the 'time_dependence' calculation type, the 
%                 function returns the observables detected using 
%                 the coil states specified at each point in time;
%                 with the 'steady_state' option specified, the 
%                 function returns the steady state values detec-
%                 ted using the coil states specified.
%
% Note: this function generates its own Liouvillian and should be 
%       called directly, without a context wrapper.
%
% ilya.kuprov@weizmann.ac.il
% alexander.karabanov@nottingham.ac.uk
% walter.kockenberger@nottingham.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=solid_effect.m>

function answer=solid_effect(spin_system,parameters)

% Check consistency
grumble(spin_system,parameters)

% Build H+ and add microwave term
report(spin_system,'solid_effect: building H+...');
[Hp,Q]=hamiltonian(assume(spin_system,'se_dnp_h+'));
Hp=Hp+orientation(Q,[0 0 0])+0.25*parameters.mw_pwr*operator(spin_system,'L-','E');

% Build H0
report(spin_system,'solid_effect: building H0...');
[H0,Q]=hamiltonian(assume(spin_system,'se_dnp_h0'));
H0=H0+orientation(Q,[0 0 0]);

% Build H- and add microwave term
report(spin_system,'solid_effect: building H-...');
[Hm,Q]=hamiltonian(assume(spin_system,'se_dnp_h-'));
Hm=Hm+orientation(Q,[0 0 0])+0.25*parameters.mw_pwr*operator(spin_system,'L+','E');

% Decide how to proceed
switch parameters.theory
    
    case 'exact'
        
        % Build the exact Hamiltonian
        H=Hp+H0+Hm+parameters.nuclear_frq*operator(spin_system,'Lz','nuclei')+...
                  -parameters.nuclear_frq*operator(spin_system,'Lz','E');
        
    otherwise
        
        % Build an average Hamiltonian
        H=average(spin_system,Hp,H0,Hm,parameters.nuclear_frq,parameters.theory);
        
end

% Set state option
if isfield(parameters,'coil')

    % Use user-specified states
    coils=parameters.coil; 
          
else
    
    % Set detection states to Lz on every spin
    coils=cell(1,spin_system.comp.nspins);
    for n=1:spin_system.comp.nspins              
        coils{n}=state(spin_system,{'Lz'},{n});        
    end
    coils=cell2mat(coils);
    
end 

% Set the initial state to thermal equilibrium
[I,Q]=hamiltonian(assume(spin_system,'labframe'),'left');
rho_eq=equilibrium(spin_system,I,Q,[0 0 0]);

% Decide the output
switch parameters.calc_type
    
    case 'time_dependence'
        
        % Get thermalized relaxation superoperator
        spin_system.rlx.equilibrium='IME';
        R=relaxation(spin_system,[0 0 0]);
        
        % Run the simulation
        answer=evolution(spin_system,H+1i*R,coils,rho_eq,parameters.time_step,...
                         parameters.n_steps,'multichannel');
                     
    case 'steady_state'
        
        % Get unthermalized relaxation superoperator
        spin_system.rlx.equilibrium='zero'; R=relaxation(spin_system,[0 0 0]);
        
        % Remove unit state singularity
        rho_unit=unit_state(spin_system); R=R-kron(rho_unit,rho_unit');
        
        % Run the simulation
        answer=evolution(spin_system,H+1i*R,coils,-R*rho_eq,parameters.time_step,...
                         parameters.n_steps,'total');
                     
    case 'trajectory' 
     
        % Get thermalized relaxation superoperator
        spin_system.rlx.equilibrium='IME'; R=relaxation(spin_system);
        
        % Run the simulation
        answer=evolution(spin_system,H+1i*R,[],rho_eq,parameters.time_step,...
                         parameters.n_steps,'trajectory');

end

end

% Consistency enforcement
function grumble(spin_system,parameters)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only available for sphten-liouv and zeeman-liouv formalisms.');
end
electron_number=find(cellfun(@(x)strcmp(x(1),'E'),spin_system.comp.isotopes));
if numel(electron_number)~=1
    error('the system must have exactly one electron.');
end
if ~isfield(parameters,'mw_pwr')
    error('microwave power should be specified in parameters.mw_pwr variable.');
end
if numel(parameters.mw_pwr)~=1
    error('parameters.mw_pwr array should have exactly one element.');
end
if ~isfield(parameters,'theory')
    error('theory level should be specified in parameters.theory variable.');
end
if ~ischar(parameters.theory)
    error('parameters.theory variable must be a character string.');
end
if ~ismember(parameters.theory,{'exact','ah_first_order','ah_second_order','ah_third_order',...
                                'kb_first_order','kb_second_order','kb_third_order','matrix_log'})
    error('incorrect parameters.theory specification, see the function header.');
end
if ~isfield(parameters,'calc_type')
    error('calculation type should be specified in parameters.calc_type variable.');
end
if ~ischar(parameters.calc_type)
    error('parameters.calc_type variable must be a character string.');
end
if ~ismember(parameters.calc_type,{'time_dependence','steady_state','trajectory'})
    error('incorrect parameters.calc_type specification, see the function header.');
end
if ~isfield(parameters,'nuclear_frq')
    error('nuclear frequency should be specified in parameters.nuclear_frq variable.');
end
if numel(parameters.nuclear_frq)~=1
    error('parameters.nuclear_frq array should have exactly one element.');
end
if strcmp(parameters.calc_type,'time_dependence')
    if ~isfield(parameters,'n_steps')
        error('number of time steps should be specified in parameters.n_steps variable.');
    end
    if numel(parameters.n_steps)~=1
        error('parameters.n_steps array should have exactly one element.');
    end
    if ~isfield(parameters,'time_step')
        error('time step length should be specified in parameters.time_step variable.');
    end
    if numel(parameters.time_step)~=1
        error('parameters.time_step array should have exactly one element.');
    end
end
end

% Mid-2011, MPLS senior staff meeting, Oxford.
%
% Tim Softley: "...and we must, of course, mention our students' role
% in this success - after all, it is they that do most of the work in
% our labs..."
%
% Incredulous voice from the audience: "Speak for yourself!"

