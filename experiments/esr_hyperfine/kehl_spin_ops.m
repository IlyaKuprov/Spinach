% defines the spin operators using Spinach operator() calls
%
% input parameters, Spinach syntax:
% spin_system: Spinach spin system object
% endor_spins: indices of ENDOR nuclei in spin_system.comp.isotopes
% N_SpinSys: number spin systems
%
% input parameters, legacy syntax:
% S: electron spin quantum number
% I: nuclear spin quantum numbers (array)
% N_Nuclei: number nuclei
% N_SpinSys: number spin systems
%
% output parameters:
% Ops: the Map with Spinach-generated spin operators
%
% February 2024 A. Kehl (akehl@gwdg.de)
% Spinach architecture migration May 2026 Talos


function Ops=kehl_spin_ops(varargin)

    % Check consistency
    grumble(varargin{:});
    if nargin>=2 && isstruct(varargin{1}) && isfield(varargin{1},'bas')
        parent_system=varargin{1};
        endor_spins=varargin{2};
        if nargin>=3
            N_SpinSys=varargin{3};
        else
            N_SpinSys=1;
        end
        N_Nuclei=numel(endor_spins);
        electron_idx=find(cellfun(@is_electron,parent_system.comp.isotopes),1);
        isotopes=parent_system.comp.isotopes([electron_idx endor_spins]);
        if N_SpinSys>1
            isotopes=isotopes(1:2);
        end
        spin_system=local_spin_system(isotopes);
    elseif nargin==1 && isa(varargin{1},'containers.Map')
        spinSys=varargin{1};
        S=spinSys("S");
        I=spinSys("I");
        N_Nuclei=spinSys("Ni_ENDOR");
        N_SpinSys=spinSys("N_SpinSys");
        if isKey(spinSys,"Nuclei")
            Nuclei=spinSys("Nuclei");
        else
            Nuclei={};
        end
        isotopes=legacy_isotopes(S,I,Nuclei,N_Nuclei,N_SpinSys);
        spin_system=local_spin_system(isotopes);
    elseif nargin==4
        S=varargin{1};
        I=varargin{2};
        N_Nuclei=varargin{3};
        N_SpinSys=varargin{4};
        isotopes=legacy_isotopes(S,I,{},N_Nuclei,N_SpinSys);
        spin_system=local_spin_system(isotopes);
    else
        error('kehl_spin_ops expects a Spinach spin system with ENDOR spin indices, a spin-system Map, or S,I,N_Nuclei,N_SpinSys.');
    end

    Ops=containers.Map;
    Ops("spin_system")=spin_system;
    Ops("isotopes")=isotopes;

    Ops("Sx")=full(operator(spin_system,'Lx',1));
    Ops("Sy")=full(operator(spin_system,'Ly',1));
    Ops("Sz")=full(operator(spin_system,'Lz',1));

    N_Nuclei_ops=numel(isotopes)-1;
    Ix=cell(1,N_Nuclei_ops);
    Iy=cell(1,N_Nuclei_ops);
    Iz=cell(1,N_Nuclei_ops);

    for n=1:N_Nuclei_ops
        spin_index=n+1;
        Ix{n}=full(operator(spin_system,'Lx',spin_index));
        Iy{n}=full(operator(spin_system,'Ly',spin_index));
        Iz{n}=full(operator(spin_system,'Lz',spin_index));
    end

    Ops("Ix")=Ix;
    Ops("Iy")=Iy;
    Ops("Iz")=Iz;
end

function spin_system=local_spin_system(isotopes)
    sys.magnet=0.35;
    sys.isotopes=isotopes;
    sys.output='hush';
    sys.disable={'hygiene'};

    inter=[];

    bas.formalism='zeeman-hilb';
    bas.approximation='none';

    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
end

function isotopes=legacy_isotopes(S,I,Nuclei,N_Nuclei,N_SpinSys)
    if N_SpinSys>1
        N_Nuclei=1;
    end

    isotopes=cell(1,N_Nuclei+1);
    isotopes{1}=electron_label(S);

    for n=1:N_Nuclei
        if numel(Nuclei)>=n
            isotope=Nuclei{n};
            if isstring(isotope)
                isotope=char(isotope);
            end
            if strcmp(isotope,'2D')
                isotope='2H';
            end
            isotopes{n+1}=isotope;
        else
            isotopes{n+1}=nucleus_label(I(n));
        end
    end
end

function label=electron_label(S)
    multiplicity=round(2*S+1);
    if abs(multiplicity-(2*S+1))>eps
        error('electron spin multiplicity must be an integer.');
    end
    if multiplicity==2
        label='E';
    else
        label=['E' num2str(multiplicity)];
    end
end

function label=nucleus_label(I)
    multiplicity=round(2*I+1);
    if multiplicity==2
        label='1H';
    elseif multiplicity==3
        label='2H';
    elseif multiplicity==6
        label='17O';
    else
        error('no default isotope label for nuclear spin multiplicity %d.',multiplicity);
    end
end

function tf=is_electron(label)
    tf=strcmp(label,'E')||~isempty(regexp(label,'^E\d+$','once'));
end

function grumble(varargin)
if nargin==0
    error('at least one input argument is required.');
end
if nargin>=2 && isstruct(varargin{1}) && isfield(varargin{1},'bas')
    if ~isnumeric(varargin{2})
        error('ENDOR spin indices must be numeric.');
    end
elseif nargin==1 && isa(varargin{1},'containers.Map')
    return
elseif nargin==4
    if (~isnumeric(varargin{1}))||(~isnumeric(varargin{2}))||(~isnumeric(varargin{3}))||(~isnumeric(varargin{4}))
        error('legacy spin-operator inputs must be numeric.');
    end
else
    error('kehl_spin_ops expects a Spinach spin system with ENDOR spin indices, a spin-system Map, or S,I,N_Nuclei,N_SpinSys.');
end
end

