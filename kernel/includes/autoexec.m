% This include is executed at the start of create.m, it over-
% rides all user input. A good use case is forcing polyadic
% or GPU arithmetic, or some other specific hardware or soft-
% ware configuration.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=autoexec.m>

% Kill stupid ass figure defaults in R2025a and later 
set(groot,'defaultFigurePosition',[680 458 560 420]); 
set(groot,'defaultFigureWindowStyle','normal'); 
set(groot,'defaultFigureMenuBar','figure'); 
set(groot,'defaultFigureToolbar','figure'); 

% IK group system settings
switch getenv('COMPUTERNAME')

    case 'ELMINSTER' % 256 AMD cores, 3 TB of RAM
        
        % 256 workers crash Matlab
        sys.parallel={'processes',128};

    case 'ALAUNDO' % 128 Intel cores, 4 TB of RAM, 8 H200 GPUs
        
        % Are GPUs involved?
        if isfield(sys,'enable')&&...
           ismember('gpu',sys.enable)

            % 4 workers per GPU are safe
            sys.parallel={'processes',32};

        else

            % Without GPUs, use all cores
            sys.parallel={'processes',128};

        end

    case 'TALOS' % 56 Intel cores, 1 TB of RAM, 3 A800 GPUs
        
        % Are GPUs involved?
        if isfield(sys,'enable')&&...
           ismember('gpu',sys.enable)

            % 4 workers per GPU are safe
            sys.parallel={'processes',12};

        else

            % Without GPUs, use all cores
            sys.parallel={'processes',56};

        end

    otherwise

        % Do nothing

end

% This relocates the scratch folder
% sys.scratch='/somewhere/it/can/write'

% ------------------------------------------------------------

% This question often produced vague answers with some candida-
% tes invoking exchange energy as if they were summoning Volde-
% mort (i.e. with fear and apprehension and/or not knowing what
% they were getting themselves into). It was also somewhat wor-
% rying to see how many of our finalists proposed a square-pla-
% nar geometry for nickel(0) complexes of the type Ni(CO)3L.
%
% Examiners' Report 2017,
% Oxford Chemistry

% #NHEAD #NGRUM