%ENDOR_KEHL_CONTEXT Spinach-style ENDOR context for the Kehl pulse sequences.
%
%   [ENDOR,ENDOR_LB,X,V]=ENDOR_KEHL_CONTEXT(SPIN_SYSTEM,PULSE_SEQUENCE,...
%   PARAMETERS,ASSUMPTIONS) prepares the Spinach-backed spin operators,
%   EPR orientation selection, ENDOR sweep axes, and line broadening, then
%   calls a pulse sequence function with Spinach's standard
%   (spin_system,parameters,H,R,K) experiment signature.
%
%   For CP ENDOR, the return signature is
%   [ENDOR,ENDOR_LB,X,Y,V].


function varargout=endor_kehl_context(spin_system,pulse_sequence,parameters,assumptions)
if nargin<4
    assumptions='labframe';
end

% Check consistency
grumble(spin_system,pulse_sequence,parameters,assumptions);


constants=kehl_constants();
spinSys=context_spin_system(spin_system,parameters,constants);
spinOps=kehl_spin_ops(spin_system,parameters.endor_spins,spinSys('N_SpinSys'));
expt=parameters.expt;
opt=parameters.opt;

if nargin>=4 && ~isempty(assumptions)
    spin_system=assume(spin_system,assumptions);
end

if isfield(parameters,'time_domain') && parameters.time_domain
    t=expt('t');
    expt('t_start')=t(6);
    paramsENDOR=kehl_prep_time(constants,spinSys,expt);
else
    paramsENDOR=kehl_prep_endor(constants,spinSys,expt);
end

paramsEPR=kehl_prep_epr(spinSys,expt);
if opt('freqDomain')==false
    epr=kehl_ori_field(constants,spinSys,spinOps,paramsEPR,paramsENDOR,opt,expt);
else
    epr=kehl_ori_freq(constants,spinSys,spinOps,paramsEPR,paramsENDOR,opt,expt);
end

localpar=parameters;
localpar.constants=constants;
localpar.spinSys=spinSys;
localpar.spinOps=spinOps;
localpar.expt=expt;
localpar.opt=opt;
localpar.paramsENDOR=paramsENDOR;
localpar.paramsEPR=paramsEPR;
localpar.epr=epr;
localpar.assumptions=assumptions;

is_cp=isfield(parameters,'cp') && parameters.cp;
if is_cp
    [y_coords,localpar]=cp_axis(localpar);
end

H=[];
R=[];
K=[];
endor_amp=pulse_sequence(spin_system,localpar,H,R,K);
endor_amp_conv=kehl_line_broaden(endor_amp,opt,localpar.expt('range_EN'));

x_coords=paramsENDOR('x_coords');
x_coords=x_coords(:,1)';
v_L=paramsENDOR('v_L');

if is_cp
    varargout={endor_amp,endor_amp_conv,x_coords,y_coords,v_L};
else
    varargout={endor_amp,endor_amp_conv,x_coords,v_L};
end
end

function spinSys=context_spin_system(spin_system,parameters,constants)
    inter=parameters.inter;
    isotopes=spin_system.comp.isotopes;
    electron_idx=find(cellfun(@is_electron,isotopes),1);
    if isempty(electron_idx)
        error('spin_system must contain an electron spin.');
    end

    endor_spins=parameters.endor_spins;
    if isfield(parameters,'epr_spins')
        epr_spins=parameters.epr_spins;
    else
        epr_spins=[];
    end
    if isfield(parameters,'n_spin_systems')
        n_spin_systems=parameters.n_spin_systems;
    else
        n_spin_systems=1;
    end

    [~,electron_mult]=spin(isotopes{electron_idx});
    n_endor=numel(endor_spins);

    spinSys=containers.Map;
    spinSys('S')=(electron_mult-1)/2;
    spinSys('g')=matrix_from_cell(inter.zeeman.matrix,electron_idx,electron_idx);
    spinSys('g_iso')=trace(spinSys('g'))/3;
    spinSys('Ni_ENDOR')=n_endor;
    spinSys('Nuclei')=isotopes(endor_spins);
    spinSys('N_SpinSys')=n_spin_systems;
    spinSys('I')=spin_numbers(isotopes(endor_spins));

    A=zeros(3*n_endor,3);
    Q=zeros(3*n_endor,3);
    Q_used=false;
    for n=1:n_endor
        spin_idx=endor_spins(n);
        A(3*n-2:3*n,:)=matrix_from_cell(inter.coupling.matrix,electron_idx,spin_idx);
        Q_block=matrix_from_cell(inter.coupling.matrix,spin_idx,spin_idx);
        Q(3*n-2:3*n,:)=Q_block;
        Q_used=Q_used||any(Q_block(:));
    end
    spinSys('A')=A;
    if isfield(parameters,'endor_quadrupole_matrix')
        Q=parameters.endor_quadrupole_matrix;
        if ~isequal(size(Q),[3*n_endor,3])
            error('parameters.endor_quadrupole_matrix must be a 3*N_ENDOR by 3 matrix.');
        end
        Q_used=any(Q(:));
    end
    spinSys('Q')=Q;
    spinSys('Q_used')=Q_used;

    CS=zeros(3*n_endor,3);
    CS_used=false;
    if isfield(inter,'zeeman') && isfield(inter.zeeman,'matrix')
        for n=1:n_endor
            spin_idx=endor_spins(n);
            if numel(inter.zeeman.matrix)>=spin_idx && ~isempty(inter.zeeman.matrix{spin_idx})
                CS_block=inter.zeeman.matrix{spin_idx}*1e6;
                CS(3*n-2:3*n,:)=CS_block;
                CS_used=CS_used||any(CS_block(:));
            end
        end
    end
    if CS_used
        spinSys('CS')=CS;
    end
    spinSys('CS_used')=CS_used;

    if isfield(parameters,'dipolar_pairs') && ~isempty(parameters.dipolar_pairs)
        pairs=parameters.dipolar_pairs;
        D=zeros(3*size(pairs,1),3);
        for n=1:size(pairs,1)
            D(3*n-2:3*n,:)=matrix_from_cell(inter.coupling.matrix,pairs(n,1),pairs(n,2));
        end
        spinSys('D')=D;
        spinSys('D_used')=true;
    else
        spinSys('D')=zeros(0,3);
        spinSys('D_used')=false;
    end

    spinSys('EPR_Nucs_used')=~isempty(epr_spins);
    if ~isempty(epr_spins)
        n_epr=numel(epr_spins);
        spinSys('EPR_Nuclei')=isotopes(epr_spins);
        spinSys('Ni_EPR')=n_epr;
        spinSys('I_EPR')=spin_numbers(isotopes(epr_spins));
        g_N_EPR=zeros(n_epr,1);
        A_EPR=zeros(3*n_epr,3);
        Q_EPR=zeros(3*n_epr,3);
        EPR_Q_used=false;
        for n=1:n_epr
            spin_idx=epr_spins(n);
            g_N_EPR(n)=kehl_nuc_gamma(constants,isotopes{spin_idx});
            A_EPR(3*n-2:3*n,:)=matrix_from_cell(inter.coupling.matrix,electron_idx,spin_idx);
            Q_block=matrix_from_cell(inter.coupling.matrix,spin_idx,spin_idx);
            Q_EPR(3*n-2:3*n,:)=Q_block;
            EPR_Q_used=EPR_Q_used||any(Q_block(:));
        end
        if isfield(parameters,'epr_quadrupole_matrix')
            Q_EPR=parameters.epr_quadrupole_matrix;
            if ~isequal(size(Q_EPR),[3*n_epr,3])
                error('parameters.epr_quadrupole_matrix must be a 3*N_EPR by 3 matrix.');
            end
            EPR_Q_used=any(Q_EPR(:));
        end
        spinSys('g_N_EPR')=g_N_EPR;
        spinSys('A_EPR')=A_EPR;
        spinSys('Q_EPR')=Q_EPR;
        spinSys('EPR_Q_used')=EPR_Q_used;
    else
        spinSys('EPR_Q_used')=false;
    end
end

function [y_coords,localpar]=cp_axis(localpar)
    expt=localpar.expt;
    opt=localpar.opt;
    if opt('powder')==false
        expt('Npts_CP')=1;
        expt('range_CP')=0;
        y_coords=1;
    else
        if expt('Npts_CP')>1
            step_CP=expt('range_CP')/(expt('Npts_CP')-1);
            y_coords=zeros(expt('Npts_CP'));
            for n=1:expt('Npts_CP')
                y_coords(n)=expt('start_CP')+(n-1)*step_CP;
            end
            y_coords=y_coords(:,1)';
        else
            y_coords=expt('start_CP');
            y_coords=y_coords(:,1)';
        end
        localpar.paramsENDOR('y_coords')=y_coords;
    end
    localpar.expt=expt;
end

function M=matrix_from_cell(cells,row,col)
    M=zeros(3,3);
    if isempty(cells)
        return
    end
    if size(cells,1)>=row && size(cells,2)>=col && ~isempty(cells{row,col})
        M=cells{row,col};
    elseif size(cells,1)>=col && size(cells,2)>=row && ~isempty(cells{col,row})
        M=cells{col,row};
    end
end

function qnums=spin_numbers(labels)
    qnums=zeros(numel(labels),1);
    for n=1:numel(labels)
        [~,multiplicity]=spin(labels{n});
        qnums(n)=(multiplicity-1)/2;
    end
end

function tf=is_electron(label)
    tf=strcmp(label,'E')||~isempty(regexp(label,'^E\d+$','once'));
end

% Consistency enforcement
function grumble(spin_system,pulse_sequence,parameters,assumptions)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isa(pulse_sequence,'function_handle')
    error('pulse_sequence must be a function handle.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if ~isfield(parameters,'inter')
    error('parameters.inter must contain the raw Spinach interaction structure.');
end
if ~isfield(parameters,'endor_spins')
    error('parameters.endor_spins must list ENDOR nuclei by Spinach spin index.');
end
if ~isfield(parameters,'expt')
    error('parameters.expt must contain the experiment parameter Map.');
end
if ~isfield(parameters,'opt')
    error('parameters.opt must contain the option Map.');
end
if (~ischar(assumptions))&&(~isstring(assumptions))
    error('assumptions must be a character string or a string scalar.');
end
end

