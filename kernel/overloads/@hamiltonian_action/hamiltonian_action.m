% Descriptor-backed matrix-free Hamiltonian action object. Syntax:
%
%        H=hamiltonian_action(spin_system,descr,euler_angles,operator_type)
%
% Parameters:
%
%    spin_system   - Spinach spin system object with a basis set and
%                    Hamiltonian assumptions
%
%    descr         - Hamiltonian descriptor table returned as the third
%                    output of hamiltonian.m
%
%    euler_angles  - a 1x3 vector specifying Euler angles, radians
%
%    operator_type - Liouville-space superoperator type: 'left',
%                    'right', 'comm', or 'acomm'
%
% Outputs:
%
%    H             - compact matrix-like object whose left multiplication
%                    by a vector returns the Hamiltonian action
%
% ilya.kuprov@weizmann.ac.il
% aditya.dev@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=hamiltonian_action.m>

classdef (InferiorClasses={?gpuArray}) hamiltonian_action

    % Object storage
    properties
        spin_system
        descr
        euler_angles
        operator_type
        dim
        coeff_iso
        coeff_aniso
        giant
        zero_tol
    end

    % Method description
    methods

        % Constructor function
        function H=hamiltonian_action(spin_system,descr,euler_angles,operator_type)

            % Check consistency
            grumble(spin_system,descr,euler_angles,operator_type);

            % Store the reduced spin system
            H.spin_system.sys=spin_system.sys;
            H.spin_system.tols=spin_system.tols;
            H.spin_system.bas=spin_system.bas;
            H.spin_system.comp=spin_system.comp;

            % Store the descriptor and control parameters
            H.descr=descr;
            H.euler_angles=euler_angles(:).';
            H.operator_type=operator_type;
            H.zero_tol=spin_system.tols.liouv_zero;

            % Get matrix dimension
            H.dim=size(mprealloc(spin_system,0),1);

            % Store the isotropic descriptor coefficients
            H.coeff_iso=descr.isotropic;

            % Preallocate anisotropic descriptor coefficients
            H.coeff_aniso=zeros(size(descr.isotropic),'like',1i);

            % Contract the spherical ranks with Wigner matrices
            for r=1:2

                % Compute the Wigner matrix
                W=wigner(r,H.euler_angles(1),...
                           H.euler_angles(2),...
                           H.euler_angles(3));

                % Add the rank contribution
                for k=1:(2*r+1)
                    for m=1:(2*r+1)
                        H.coeff_aniso=H.coeff_aniso+W(k,m)*...
                                       descr.ist_coeff(:,r^2+k-1).*...
                                       descr.irr_comp(:,r^2+m-1);
                    end
                end

            end

            % Preallocate the oriented giant spin contribution
            H.giant=spalloc(H.dim,H.dim,0);

            % Build the oriented giant spin contribution
            for n=1:spin_system.comp.nspins
                for r=1:numel(spin_system.inter.giant.coeff{n})

                    % Compute the Wigner matrix
                    W=wigner(r,H.euler_angles(1),...
                               H.euler_angles(2),...
                               H.euler_angles(3));

                    % Process the giant spin assumption
                    switch spin_system.inter.giant.strength{n}

                        case 'strong'

                            % Loop over spherical tensor indices
                            for k=1:(2*r+1)

                                % Contract coefficient components
                                coeff=W(k,:)*spin_system.inter.giant.coeff{n}{r}(:);

                                % Add the operator contribution
                                if abs(coeff)>H.zero_tol
                                    ist_spec=['T' num2str(r) ',' num2str(r-k+1)];
                                    xyz=operator(H.spin_system,{ist_spec},{n},...
                                                 H.operator_type,'xyz');
                                    H.giant=H.giant+sparse(xyz(:,1),xyz(:,2),...
                                                           coeff*xyz(:,3),H.dim,H.dim);
                                end

                            end

                        case 'secular'

                            % Contract coefficient components
                            coeff=W(r+1,:)*spin_system.inter.giant.coeff{n}{r}(:);

                            % Add the operator contribution
                            if abs(coeff)>H.zero_tol
                                ist_spec=['T' num2str(r) ',0'];
                                xyz=operator(H.spin_system,{ist_spec},{n},...
                                             H.operator_type,'xyz');
                                H.giant=H.giant+sparse(xyz(:,1),xyz(:,2),...
                                                       coeff*xyz(:,3),H.dim,H.dim);
                            end

                        case 'ignore'

                            % Move to the next term
                            continue

                        otherwise

                            % Bomb out with unexpected strength parameters
                            error('unknown giant spin interaction strength specification.');

                    end

                end
            end

            % Match orientation.m Hermitisation
            H.giant=(H.giant+H.giant')/2;

        end

        % Numeric property
        function answer=isnumeric(H) %#ok<MANU>

            answer=true();

        end

        % Matrix property
        function answer=ismatrix(H) %#ok<MANU>

            answer=true();

        end

        % Number of elements
        function answer=numel(H)

            answer=H.dim^2;

        end

        % Usually not unit matrix
        function answer=iseye(H) %#ok<MANU>

            answer=false();

        end

    end

end

% Consistency enforcement
function grumble(spin_system,descr,euler_angles,operator_type)
if ~isfield(spin_system,'bas')
    error('basis set information is missing, run basis() before calling this function.');
end
if ~ismember(spin_system.bas.formalism,{'zeeman-liouv','sphten-liouv'})
    error('hamiltonian_action only supports Liouville-space formalisms.');
end
if ~istable(descr)
    error('descr must be a Hamiltonian descriptor table.');
end
if ~all(ismember({'nL','nS','opL','opS','isotropic','ist_coeff','irr_comp'},...
                 descr.Properties.VariableNames))
    error('descriptor table is missing Hamiltonian descriptor fields.');
end
if (~isnumeric(descr.nL))||(~isnumeric(descr.nS))||...
   any(mod(descr.nL,1)~=0)||any(mod(descr.nS,1)~=0)||...
   any(descr.nL<0)||any(descr.nS<0)
    error('descriptor spin numbers must be non-negative integers.');
end
if (~iscell(descr.opL))||(~iscell(descr.opS))||...
   any(~cellfun(@ischar,descr.opL))||any(~cellfun(@ischar,descr.opS))
    error('descriptor operator labels must be character strings.');
end
if (~isnumeric(descr.isotropic))||(~iscolumn(descr.isotropic))
    error('descriptor isotropic coefficients must be a numeric column.');
end
if (~isnumeric(descr.ist_coeff))||(size(descr.ist_coeff,2)~=8)
    error('descriptor IST coefficients must have eight columns.');
end
if (~isnumeric(descr.irr_comp))||(size(descr.irr_comp,2)~=8)
    error('descriptor irreducible coefficients must have eight columns.');
end
if (~isnumeric(euler_angles))||(~isreal(euler_angles))||(numel(euler_angles)~=3)
    error('euler_angles must be a three-element real vector.');
end
if (~ischar(operator_type))||...
   (~ismember(operator_type,{'left','right','comm','acomm'}))
    error('incorrect operator type specification.');
end
end

