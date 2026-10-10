% Microwave frequency scan steady-state DNP experiment. Returns the
% steady-state population of the user-specified states as a function
% of microwave irradiation frequency. Syntax:
%
%          dnp=dnp_freq_scan(spin_system,parameters,H,R,K)
%
% Parameters:
%
%    parameters.mw_pwr       -   microwave power, rad/s
%
%    parameters.mw_frq       -   row vector of microwave frequ-
%                                ency offsets (rad/s) relative
%                                to the reference g-factor
%
%    parameters.g_ref        -   reference g-factor around which
%                                frequency offsets are specified
%
%    parameters.rho0         -   thermal equilibrium state
%
%    parameters.coil         -   coil state vector or a horizon-
%                                tal stack thereof
%
%    parameters.mw_oper      -   microwave irradiation operator
%
%    parameters.ez_oper      -   Lz operator on the electrons
%
%    parameters.method       -   calculation method: 'fp-backs',
%                                'fp-gmres', 'lvn-backs', or
%                                'lvn-gmres'
%
%    parameters.nphases      -   number of microwave phase grid
%                                points for the Fokker-Planck path
%
%    H - Hamiltonian matrix, received from context function
%
%    R - relaxation superoperator, received from context function
%
%    K - kinetics superoperator, received from context function
%
% Note: segmented substance inputs are not yet supported.
%
% Outputs:
%
%    dnp    -  an array of steady state expectation values for
%              the states specified in parameters.coil at each
%              of the microwave frequencies supplied
%
% Note: the relaxation superoperator must NOT be thermalized for
%       this type of calculation (inter.equilibrium='zero').
%
% With polyadics enabled, fp-gmres uses implicit FFT differentiation
% and a static Fourier-block preconditioner. Direct backslash and
% LvN methods keep their explicit implementations.
%
% ilya.kuprov@weizmann.ac.il
% alexander.karabanov@nottingham.ac.uk
% walter.kockenberger@nottingham.ac.uk
% mariagrazia.concilio@sjtu.edu.cn
%
% <https://spindynamics.org/wiki/index.php?title=dnp_freq_scan.m>

function dnp=dnp_freq_scan(spin_system,parameters,H,R,K) 

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Damp the trace direction and check
switch spin_system.bas.formalism

    case 'sphten-liouv'

        % Unit state carries the trace
        R(1,1)=-mean(abs(diag(R)));

    case 'zeeman-liouv'

        % Trace direction is the stretched unit matrix
        dim=sqrt(size(R,1)); u0=speye(dim); u0=u0(:);
        R=R-(mean(abs(diag(R)))/dim)*(u0*u0');

end
Rc=condest(R);
if Rc>1e9, error('R must be non-singular.'); end
if any(ismember({'redfield','naka-zwan'},spin_system.rlx.theories))&&...
           (Rc>1/spin_system.tols.rlx_integration)
    error('R acuracy too low, reduce spin_system.tols.rlx_integration');
end

% Preallocate the answer
dnp=zeros([numel(parameters.mw_frq) size(parameters.coil,2)],'like',1i);

switch parameters.method

    % Fokker-Planck path    
    case {'fp-backs','fp-gmres'}

        % Translate frequency offsets
        carrier_freq=-spin('E')*(parameters.g_ref/spin_system.tols.freeg)*spin_system.inter.magnet;
        parameters.mw_frq=carrier_freq+parameters.mw_frq;

        % Get problem dimensions
        spc_dim=parameters.nphases; spn_dim=size(H,1);

        % Select the implicit derivative for the iterative solver only
        use_poly=ismember('polyadic',spin_system.sys.enable)&&strcmp(parameters.method,'fp-gmres');
        deriv_system=spin_system;
        if ~use_poly
            deriv_system.sys.enable=setdiff(spin_system.sys.enable,{'polyadic'});
        end
        [phases,d_dphi]=fourdif(deriv_system,spc_dim,1);

        % Keep the spin-sized static generator for Fourier-block solves
        static_gen=H+1i*R+1i*K;
        wave_nums=fftdiff(1,spc_dim,2*pi/spc_dim);
        if mod(spc_dim,2)==0, wave_nums(spc_dim/2+1)=0; end

        % Lift operators without opening implicit Fourier factors
        cosines=spdiags(cos(phases),0,spc_dim,spc_dim);
        if use_poly
            PT=polyadic({{d_dphi,opium(spn_dim,1)}});
            H=polyadic({{opium(spc_dim,1),H}});
            R=polyadic({{opium(spc_dim,1),R}});
            K=polyadic({{opium(spc_dim,1),K}});
            MW=polyadic({{cosines,parameters.mw_pwr*parameters.mw_oper}});
        else
            PT=kron(d_dphi,speye(spn_dim,spn_dim));
            H=kron(speye(spc_dim,spc_dim),H);
            R=kron(speye(spc_dim,spc_dim),R);
            K=kron(speye(spc_dim,spc_dim),K);
            MW=parameters.mw_pwr*kron(cosines,parameters.mw_oper);
        end

        % Precompute the right hand side
        if use_poly
            b=R*full(kron(ones(spc_dim,1),parameters.rho0));
        else
            b=R*kron(ones(spc_dim,1),parameters.rho0);
        end

        % MDCS diagnostics
        parallel_profiler_start;

        % Loop in parallel over the frequencies
        parfor n=1:numel(parameters.mw_frq)

            % Build the evolution generator without materialising FFT cores
            L=H+1i*R+1i*K+MW+1i*parameters.mw_frq(n)*PT; %#ok<PFBNS>
            if ~use_poly, L=sparse(L); end

            % Get the steady state DNP
            switch parameters.method

                case 'fp-gmres'

                    % Precondition implicit GMRES with spin-sized Fourier blocks
                    if use_poly
                        inverses=cell(1,spc_dim);
                        for phase=1:spc_dim
                            inverses{phase}=decomposition(-1i*static_gen+...
                                parameters.mw_frq(n)*wave_nums(phase)*speye(spn_dim),'lu'); %#ok<PFBNS>
                        end
                        [big_inv,flag,relres]=gmres(@(x)-1i*(L*x),b,min(10,numel(b)),...
                            1e-10,numel(b),@(x)phase_inverse(x,inverses,spn_dim,spc_dim));
                        if flag~=0
                            error('implicit DNP GMRES failed: flag=%d, relative residual=%g.',flag,relres);
                        end
                    else
                        big_inv=gmres(-1i*L,b,10,1e-10,numel(b),L);
                    end

                case 'fp-backs'

                    % Use Matlab's backslash
                    big_inv=-(1i*L)\b;

                otherwise

                    % Complain and bomb out
                    big_inv=0; error('unknown solver.'); %#ok<NASGU>

            end

            % Fold phase dimension and compute observables
            dnp(n,:)=parameters.coil'*mean(reshape(big_inv,[spn_dim spc_dim]),2);

        end

        % Get MDCS diagnostics
        parallel_profiler_report;

    % Liouville space path
    case {'lvn-gmres','lvn-backs'}

        % Translate frequency offsets
        parameters.mw_frq=parameters.mw_frq-spin('E')*spin_system.inter.magnet*...
                         (parameters.g_ref-spin_system.tols.freeg)/spin_system.tols.freeg;

        % Compose the Liouvillian 
        L=sparse(H+1i*R+1i*K+(parameters.mw_pwr/2)*parameters.mw_oper);

        % Precompute the right hand side
        b=R*parameters.rho0;

        % MDCS diagnostics
        parallel_profiler_start;

        % Loop in parallel over the frequencies
        parfor n=1:numel(parameters.mw_frq)

            % Add electron offset terms to the Hamiltonian
            L_curr=sparse(L-parameters.mw_frq(n)*parameters.ez_oper); %#ok<PFBNS>

            % Get the steady state DNP
            switch parameters.method

                case 'lvn-gmres'

                    % Get the preconditioner
                    [M1,M2]=ilu(L_curr,struct('type','crout','droptol',1e-6));

                    % Run using preconditioned GMRES
                    dnp(n,:)=parameters.coil'*gmres(-1i*L_curr,b,10,1e-10,numel(b),M1,M2);

                case 'lvn-backs'

                    % Run using Matlab's backslash
                    dnp(n,:)=parameters.coil'*(-1i*L_curr\b);

                otherwise

                    % Complain and bomb out
                    error('unknown solver.');

            end

        end

        % Get MDCS diagnostics
        parallel_profiler_report;

    otherwise

        % Compain and bomb out
        error('unknown calculation method.');

end

end

% Apply the inverse static generator in the Fourier phase basis
function answer=phase_inverse(rhs,inverses,spn_dim,spc_dim)

% Transform the phase coordinate while keeping spin states contiguous
answer=fft(reshape(rhs,spn_dim,spc_dim),[],2);

% Solve the independent spin-sized Fourier blocks
for phase=1:spc_dim
    answer(:,phase)=inverses{phase}\answer(:,phase);
end

% Return to the physical phase coordinate and vectorise
answer=ifft(answer,[],2); answer=answer(:);

end

% Consistency checking
function grumble(spin_system,parameters,H,R,K)
if spin_system.bas.nsubst>1
    error('Spinach:dnp_freq_scan:segmentedSubstances',...
          'DNP scans require a single substance.');
end
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only available for sphten-liouv and zeeman-liouv formalisms.');
end
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~isfield(parameters,'mw_pwr')
    error('microwave power should be specified in parameters.mw_pwr');
elseif numel(parameters.mw_pwr)~=1
    error('parameters.mw_pwr array should have exactly one element.');
end
if ~isfield(parameters,'mw_frq')
    error('microwave frequencies should be specified in parameters.mw_frq');
elseif (~isnumeric(parameters.mw_frq))||(~isreal(parameters.mw_frq))
    error('elements of parameters.mw_frq variable must be real numbers.');
end
if ~isfield(parameters,'g_ref')
    error('reference g-factor should be specified in parameters.g_ref');
elseif (~isnumeric(parameters.g_ref))||(~isreal(parameters.g_ref))||(~isscalar(parameters.g_ref))
    error('parameters.g_ref should be a real scalar.');
end
if ~isfield(parameters,'rho0')
    error('thermal equilibrium state must be specified in parameters.rho0');
end
if ~isfield(parameters,'coil')
    error('one or more detection states must be specified in parameters.coil');
end
if ~isfield(parameters,'mw_oper')      
    error('microwave irradiation operator must be specified in parameters.mw_oper');   
end
if ~isfield(parameters,'method')
    error('calculation method must be specified in parameters.method');
end
if ~ischar(parameters.method)
    error('parameters.method must be a character string.');
end
if ismember(parameters.method,{'lvn-gmres','lvn-backs'})  
    if ~isfield(parameters,'ez_oper')
        error('electron Lz operator must be provided in parameters.ez_oper');
    end
    if ~strcmp(spin_system.inter.assumptions,'esr')
        error('LvN formalism requires esr assumptions.');
    end
end
if ismember(parameters.method,{'fp-backs','fp-gmres'})  
    if ~isfield(parameters,'nphases')
        error('MW phase grid size must be specified in parameters.nphases');
    end
    if (~isnumeric(parameters.nphases))||(~isreal(parameters.nphases))||...
       (~isscalar(parameters.nphases))||(mod(parameters.nphases,1)~=0)||...
       (parameters.nphases<1)
        error('parameters.nphases must be a positive real integer.');
    end
    if ~strcmp(spin_system.inter.assumptions,'labframe')
        error('FP formalism requires labframe assumptions.');
    end
end
end

% Too many people want to *have written*.
%
% Terry Pratchett


