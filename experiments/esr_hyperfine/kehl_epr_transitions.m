% Liouville-space EPR transition selection for Kehl ENDOR. Syntax:
%
%      [positions,moments]=kehl_epr_transitions(spin_system,...
%                           parameters,euler_angles,mode,probe)
%
% Parameters:
%
%   spin_system      - Zeeman-Liouville Spinach spin system.
%   parameters       - Kehl ENDOR context parameter structure.
%   euler_angles     - Kehl orientation angles stored by the context.
%   mode             - 'field' for transition fields, or 'frequency' for transition frequencies.
%   probe            - microwave frequency in Hz for 'field', or magnetic field in T for 'frequency'.
%
% Outputs:
%
%   positions        - transition fields in T, or transition frequencies in Hz.
%   moments          - microwave transition moments.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_epr_transitions.m>

function [positions,moments]=kehl_epr_transitions(spin_system,...
        parameters,euler_angles,mode,probe)

    % Check consistency
    grumble(spin_system,parameters,euler_angles,mode,probe);

    % Assemble field-dependent and field-independent Liouvillians
    [Hz,Hc,Hmw]=epr_liouvillian(spin_system,parameters,euler_angles);

    % Solve the transition-selection problem in Liouville space
    switch mode
        case 'field'
            omega=2*pi*probe;
            dim=size(Hc,1);
            [uv,roots]=eig(omega*speye(dim)-full(Hc),full(Hz),'vector');
            hit_list=(~isfinite(roots))|...
                (abs(imag(roots))>sqrt(eps)*max(1,abs(real(roots))));
            positions=real(roots);
            hit_list=hit_list|(positions<=0);
            if isfield(parameters,'paramsEPR')&&...
                    isKey(parameters.paramsEPR,'fieldAxis')
                field_axis=parameters.paramsEPR('fieldAxis');
                hit_list=hit_list|(positions<min(field_axis))|...
                    (positions>max(field_axis));
            end
        case 'frequency'
            H=probe*Hz+Hc;
            [uv,omega]=eig(full(H),'vector');
            hit_list=(~isfinite(omega))|...
                (abs(imag(omega))>sqrt(eps)*max(1,abs(real(omega))));
            positions=real(omega)/(2*pi);
            hit_list=hit_list|(positions<=0);
        otherwise
            error('unexpected transition-selection mode.');
    end

    % Prune unphysical roots before moment calculation
    positions(hit_list)=[];
    uv(:,hit_list)=[];

    % Return immediately if no transitions survive
    if isempty(positions)
        moments=[];
        return
    end

    % Get root-merging tolerance after invalid roots are gone
    if strcmp(mode,'field')
        tol=field_tol(parameters,positions);
    else
        tol=freq_tol(parameters,positions);
    end

    % Normalise transition dyadics
    norms=sqrt(sum(abs(uv).^2,1));
    hit_list=(norms==0)|(~isfinite(norms));
    positions(hit_list)=[];
    uv(:,hit_list)=[];
    norms(hit_list)=[];
    if isempty(positions)
        moments=[];
        return
    end
    uv=uv./norms;

    % Compute microwave transition moments
    moments=abs(Hmw'*uv).^2;

    % Prune numerically dark transitions
    hit_list=(~isfinite(moments))|...
        (moments<=max(moments)*1e-12);
    positions(hit_list)=[];
    moments(hit_list)=[];

    % Merge numerically degenerate transition roots
    [positions,moments]=merge_roots(positions,moments,tol);

end

% Assemble Zeeman and coupling Liouvillians for EPR selection
function [Hz,Hc,Hmw]=epr_liouvillian(spin_system,parameters,euler_angles)

    % Convert stored Kehl angles into the legacy direction matrix
    theta=-euler_angles(2);
    phi=-euler_angles(3);
    R=zeros(3,3);
    R(1,1)=cos(theta)*cos(phi);
    R(1,2)=cos(theta)*sin(phi);
    R(1,3)=-sin(theta);
    R(2,1)=-sin(phi);
    R(2,2)=cos(phi);
    R(3,1)=sin(theta)*cos(phi);
    R(3,2)=sin(theta)*sin(phi);
    R(3,3)=cos(theta);

    % Initialise Liouville-space Hamiltonian components
    dim=prod(spin_system.comp.mults)^2;
    Hz=sparse(dim,dim);
    Hc=sparse(dim,dim);

    % Add the orientation-projected electron Zeeman term per Tesla
    electron_idx=parameters.electron_spin_idx;
    g_lab=R*parameters.g_matrix*R';
    Sz=operator(spin_system,'Lz',electron_idx);
    Hz=Hz+2*pi*parameters.constants('MU_B')*g_lab(3,3)*...
        Sz/parameters.constants('H');

    % Add nuclear Zeeman, hyperfine, and quadrupolar terms
    for n=1:parameters.n_epr
        spin_idx=parameters.epr_spins(n);
        Iz=operator(spin_system,'Lz',spin_idx);
        Hz=Hz-2*pi*parameters.epr_gamma_hz_t(n)*Iz;
        A=R*parameters.epr_hfc_matrix(3*n-2:3*n,:)*R';
        Hc=Hc+2*pi*A(3,3)*kehl_product_comm(spin_system,...
            {'Lz','Lz'},[electron_idx spin_idx]);
        Hc=Hc+2*pi*sqrt(A(1,3)^2+A(2,3)^2)*...
            kehl_product_comm(spin_system,{'Lz','Lx'},...
            [electron_idx spin_idx]);
        if parameters.epr_nqi_active
            Q=R*parameters.epr_nqi_matrix(3*n-2:3*n,:)*R';
            Hc=Hc+2*pi*Q(1,1)*kehl_product_comm(spin_system,...
                {'Lx','Lx'},[spin_idx spin_idx]);
            Hc=Hc+2*pi*Q(1,2)*kehl_product_comm(spin_system,...
                {'Lx','Ly'},[spin_idx spin_idx]);
            Hc=Hc+2*pi*Q(1,3)*kehl_product_comm(spin_system,...
                {'Lx','Lz'},[spin_idx spin_idx]);
            Hc=Hc+2*pi*Q(2,1)*kehl_product_comm(spin_system,...
                {'Ly','Lx'},[spin_idx spin_idx]);
            Hc=Hc+2*pi*Q(2,2)*kehl_product_comm(spin_system,...
                {'Ly','Ly'},[spin_idx spin_idx]);
            Hc=Hc+2*pi*Q(2,3)*kehl_product_comm(spin_system,...
                {'Ly','Lz'},[spin_idx spin_idx]);
            Hc=Hc+2*pi*Q(3,1)*kehl_product_comm(spin_system,...
                {'Lz','Lx'},[spin_idx spin_idx]);
            Hc=Hc+2*pi*Q(3,2)*kehl_product_comm(spin_system,...
                {'Lz','Ly'},[spin_idx spin_idx]);
            Hc=Hc+2*pi*Q(3,3)*kehl_product_comm(spin_system,...
                {'Lz','Lz'},[spin_idx spin_idx]);
        end
    end

    % Build microwave transition vector
    Hmw=state(spin_system,'Lx',electron_idx);

end

% Commutation superoperator for an ordered product of spin operators
function A=kehl_product_comm(spin_system,operators,spins)

    % Initialise left and right product superoperators
    dim=prod(spin_system.comp.mults)^2;
    left_op=speye(dim);
    right_op=speye(dim);

    % Build the left product in physical operator order
    for n=1:numel(operators)
        left_op=left_op*operator(spin_system,operators{n},spins(n),'left');
    end

    % Build the right product in reverse order
    for n=numel(operators):-1:1
        right_op=right_op*operator(spin_system,operators{n},spins(n),'right');
    end

    % Return the commutation superoperator
    A=left_op-right_op;

end

% Merge transition roots that are identical to numerical precision
function [positions,moments]=merge_roots(positions,moments,tol)

    % Sort roots before merging
    [positions,idx]=sort(positions(:));
    moments=moments(:);
    moments=moments(idx);

    % Walk the sorted root list
    pos_out=[];
    mom_out=[];
    while ~isempty(positions)
        hit_list=abs(positions-positions(1))<=tol;
        root_pos=positions(hit_list);
        root_mom=moments(hit_list);
        if sum(root_mom)>0
            pos_out(end+1)=sum(root_pos.*root_mom)/sum(root_mom); %#ok<AGROW>
        else
            pos_out(end+1)=mean(root_pos); %#ok<AGROW>
        end
        mom_out(end+1)=sum(root_mom); %#ok<AGROW>
        positions(hit_list)=[];
        moments(hit_list)=[];
    end

    % Return row vectors
    positions=pos_out;
    moments=mom_out;

end

% Frequency root tolerance
function tol=freq_tol(parameters,positions)
    if isfield(parameters,'epr_freq_step_hz')
        tol=max(parameters.epr_freq_step_hz*1e-6,...
            100*eps(max(abs(positions))));
    else
        tol=max(1e-6,100*eps(max(abs(positions))));
    end
end

% Field root tolerance
function tol=field_tol(parameters,positions)
    if isfield(parameters,'field_step_t')
        tol=max(parameters.field_step_t*1e-6,...
            100*eps(max(abs(positions))));
    else
        tol=max(1e-12,100*eps(max(abs(positions))));
    end
end

% Consistency enforcement
function grumble(spin_system,parameters,euler_angles,mode,probe)
    if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
            (~isfield(spin_system,'comp'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if ~strcmp(spin_system.bas.formalism,'zeeman-liouv')
        error('spin_system must use zeeman-liouv formalism.');
    end
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
    if (~isnumeric(euler_angles))||(~isreal(euler_angles))||...
            (numel(euler_angles)~=3)
        error('euler_angles must be a three-element real vector.');
    end
    if (~ischar(mode))||(~ismember(mode,{'field','frequency'}))
        error('mode must be either ''field'' or ''frequency''.');
    end
    if (~isnumeric(probe))||(~isreal(probe))||(~isscalar(probe))
        error('probe must be a real numeric scalar.');
    end
end

