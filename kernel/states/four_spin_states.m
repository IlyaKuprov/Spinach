% Returns user-specified states for a system of four spin-1/2
% particles; see also the enclosed Mathematica file. Syntax:
%
%     rho=four_spin_states(spin_system,spins,spin_state)
%
% Arguments:
%
%   spins      - a row vector of four spin numbers
%
%   spin_state - one of the possible singlet-triplet
%                product states, see below
%
% Outputs:
%
%   rho     - a density matrix (Hilbert space) or
%             a state vector (Liouville space) 
%
% ilya.kuprov@weizmann.ac.il
% ttheis@ncsu.edu
%
% <https://spindynamics.org/wiki/index.php?title=four_spin_states.m>

function rho=four_spin_states(spin_system,spins,spin_state)

% Check consistency
grumble(spin_system,spins,spin_state);

% Component operators: four-spin
XXXX=state(spin_system,{'Lx','Lx','Lx','Lx'},num2cell(spins(:)));
XXYY=state(spin_system,{'Lx','Lx','Ly','Ly'},num2cell(spins(:)));
XXZZ=state(spin_system,{'Lx','Lx','Lz','Lz'},num2cell(spins(:)));
YYXX=state(spin_system,{'Ly','Ly','Lx','Lx'},num2cell(spins(:)));
YYYY=state(spin_system,{'Ly','Ly','Ly','Ly'},num2cell(spins(:)));
YYZZ=state(spin_system,{'Ly','Ly','Lz','Lz'},num2cell(spins(:)));
ZZXX=state(spin_system,{'Lz','Lz','Lx','Lx'},num2cell(spins(:)));
ZZYY=state(spin_system,{'Lz','Lz','Ly','Ly'},num2cell(spins(:)));
ZZZZ=state(spin_system,{'Lz','Lz','Lz','Lz'},num2cell(spins(:)));

% Component operators: three-spin
XXEZ=state(spin_system,{'Lx','Lx','E' ,'Lz'},num2cell(spins(:)));
XXZE=state(spin_system,{'Lx','Lx','Lz','E' },num2cell(spins(:)));
YYEZ=state(spin_system,{'Ly','Ly','E' ,'Lz'},num2cell(spins(:)));
YYZE=state(spin_system,{'Ly','Ly','Lz','E' },num2cell(spins(:)));
ZZEZ=state(spin_system,{'Lz','Lz','E' ,'Lz'},num2cell(spins(:)));
ZZZE=state(spin_system,{'Lz','Lz','Lz','E' },num2cell(spins(:)));
EZXX=state(spin_system,{'E' ,'Lz','Lx','Lx'},num2cell(spins(:)));
EZYY=state(spin_system,{'E' ,'Lz','Ly','Ly'},num2cell(spins(:)));
EZZZ=state(spin_system,{'E' ,'Lz','Lz','Lz'},num2cell(spins(:)));
ZEXX=state(spin_system,{'Lz','E' ,'Lx','Lx'},num2cell(spins(:)));
ZEYY=state(spin_system,{'Lz','E' ,'Ly','Ly'},num2cell(spins(:)));
ZEZZ=state(spin_system,{'Lz','E' ,'Lz','Lz'},num2cell(spins(:)));

% Component operators: two-spin
EEXX=state(spin_system,{'E', 'E' ,'Lx','Lx'},num2cell(spins(:)));
EEYY=state(spin_system,{'E', 'E' ,'Ly','Ly'},num2cell(spins(:)));
EEZZ=state(spin_system,{'E', 'E' ,'Lz','Lz'},num2cell(spins(:)));
XXEE=state(spin_system,{'Lx','Lx','E' ,'E' },num2cell(spins(:)));
YYEE=state(spin_system,{'Ly','Ly','E' ,'E' },num2cell(spins(:)));
ZZEE=state(spin_system,{'Lz','Lz','E' ,'E' },num2cell(spins(:)));
EZEZ=state(spin_system,{'E' ,'Lz','E' ,'Lz'},num2cell(spins(:)));
EZZE=state(spin_system,{'E' ,'Lz','Lz','E' },num2cell(spins(:)));
ZEEZ=state(spin_system,{'Lz','E', 'E' ,'Lz'},num2cell(spins(:)));
ZEZE=state(spin_system,{'Lz','E', 'Lz','E' },num2cell(spins(:)));

% Component operators: one-spin
EEEZ=state(spin_system,{'E', 'E' ,'E' ,'Lz'},num2cell(spins(:)));
EEZE=state(spin_system,{'E', 'E' ,'Lz','E' },num2cell(spins(:)));
EZEE=state(spin_system,{'E', 'Lz','E' ,'E' },num2cell(spins(:)));
ZEEE=state(spin_system,{'Lz','E' ,'E' ,'E' },num2cell(spins(:)));

% Component operators: unit state
EEEE=state(spin_system,{'E', 'E', 'E', 'E' },num2cell(spins(:)));

% Build the state
switch spin_state

    case 'S(x)S'   % Singlet (x) singlet
        
        rho=EEEE/16-(EEXX+EEYY+EEZZ+XXEE+YYEE+ZZEE)/4+...
            XXXX+XXYY+XXZZ+YYXX+YYYY+YYZZ+ZZXX+ZZYY+ZZZZ;

    case 'S(x)TU'  % Singlet (x) triplet-up

        rho=EEEE/16+(EEEZ+EEZE)/8+(EEZZ-XXEE-YYEE-ZZEE)/4-...
            (XXEZ+XXZE+YYEZ+YYZE+ZZEZ+ZZZE)/2-XXZZ-YYZZ-ZZZZ;

    case 'S(x)T0'  % Singlet (x) triplet-middle

        rho=EEEE/16+(EEXX+EEYY-EEZZ-XXEE-YYEE-ZZEE)/4-...
            XXXX-XXYY+XXZZ-YYXX-YYYY+YYZZ-ZZXX-ZZYY+ZZZZ;
        
    case 'S(x)TD'  % Singlet (x) triplet-down
        
        rho=EEEE/16-(EEEZ+EEZE)/8+(EEZZ-XXEE-YYEE-ZZEE)/4+...
            (XXEZ+XXZE+YYEZ+YYZE+ZZEZ+ZZZE)/2-XXZZ-YYZZ-ZZZZ;

    case 'TU(x)S'  % Triplet-up (x) singlet
        
        rho=EEEE/16+(EZEE+ZEEE)/8-(EEXX+EEYY+EEZZ-ZZEE)/4-...
            (EZXX+EZYY+EZZZ+ZEXX+ZEYY+ZEZZ)/2-ZZXX-ZZYY-ZZZZ;

    case 'T0(x)S'  % Triplet-middle (x) singlet
        
        rho=EEEE/16-(EEXX+EEYY+EEZZ-XXEE-YYEE+ZZEE)/4-...
            XXXX-XXYY-XXZZ-YYXX-YYYY-YYZZ+ZZXX+ZZYY+ZZZZ;
        
    case 'TD(x)S'  % Triplet-down (x) singlet
        
        rho=EEEE/16-(EZEE+ZEEE)/8-(EEXX+EEYY+EEZZ-ZZEE)/4+...
            (EZXX+EZYY+EZZZ+ZEXX+ZEYY+ZEZZ)/2-ZZXX-ZZYY-ZZZZ;

    case 'TU(x)TU' % Triplet-up (x) triplet-up
        
        rho=EEEE/16+(EEEZ+EEZE+EZEE+ZEEE)/8+...
            (EEZZ+EZEZ+EZZE+ZEEZ+ZEZE+ZZEE)/4+...
            (EZZZ+ZEZZ+ZZEZ+ZZZE)/2+ZZZZ;

    case 'T0(x)TU' % Triplet-middle (x) triplet-up

        rho=EEEE/16+(EEEZ+EEZE)/8+(EEZZ+XXEE+YYEE-ZZEE)/4+...
            (XXEZ+XXZE+YYEZ+YYZE-ZZEZ-ZZZE)/2+XXZZ+YYZZ-ZZZZ;
        
    case 'TD(x)TU' % Triplet-down (x) triplet-up
        
        rho=EEEE/16+(EEEZ+EEZE-EZEE-ZEEE)/8+...
            (EEZZ-EZEZ-EZZE-ZEEZ-ZEZE+ZZEE)/4-...
            (EZZZ+ZEZZ-ZZEZ-ZZZE)/2+ZZZZ;

    case 'TU(x)T0' % Triplet-up (x) triplet-middle
        
        rho=EEEE/16+(EZEE+ZEEE)/8+(EEXX+EEYY-EEZZ+ZZEE)/4+...
            (EZXX+EZYY-EZZZ+ZEXX+ZEYY-ZEZZ)/2+ZZXX+ZZYY-ZZZZ;

    case 'T0(x)T0' % Triplet-middle (x) triplet-middle
        
        rho=EEEE/16+(EEXX+EEYY-EEZZ+XXEE+YYEE-ZZEE)/4+...
            XXXX+XXYY-XXZZ+YYXX+YYYY-YYZZ-ZZXX-ZZYY+ZZZZ;
        
    case 'TD(x)T0' % Triplet-down (x) triplet-middle
        
        rho=EEEE/16-(EZEE+ZEEE)/8+(EEXX+EEYY-EEZZ+ZZEE)/4-...
            (EZXX+EZYY-EZZZ+ZEXX+ZEYY-ZEZZ)/2+ZZXX+ZZYY-ZZZZ;

    case 'TU(x)TD' % Triplet-up (x) triplet-down
        
        rho=EEEE/16-(EEEZ+EEZE-EZEE-ZEEE)/8+...
            (EEZZ-EZEZ-EZZE-ZEEZ-ZEZE+ZZEE)/4+...
            (EZZZ+ZEZZ-ZZEZ-ZZZE)/2+ZZZZ;

    case 'T0(x)TD' % Triplet-middle (x) triplet-down
        
        rho=EEEE/16-(EEEZ+EEZE)/8+(EEZZ+XXEE+YYEE-ZZEE)/4-...
            (XXEZ+XXZE+YYEZ+YYZE-ZZEZ-ZZZE)/2+XXZZ+YYZZ-ZZZZ;
        
    case 'TD(x)TD' % Triplet-down (x) triplet-down
        
        rho=EEEE/16-(EEEZ+EEZE+EZEE+ZEEE)/8+...
            (EEZZ+EZEZ+EZZE+ZEEZ+ZEZE+ZZEE)/4-...
            (EZZZ+ZEZZ+ZZEZ+ZZZE)/2+ZZZZ;

    otherwise

        % Complain and bomb out
        error('unknown spin state.');

end

end

% Consistency enforcement
function grumble(spin_system,spins,spin_state)
if (~isnumeric(spins))||(~isvector(spins))||(numel(spins)~=4)||...
   (~isreal(spins))||any(spins<1,'all')||any(mod(spins,1)~=0,'all')||...
   (numel(unique(spins))~=numel(spins))||(~isrow(spins))
    error('spins must be a row vector of four different real integers.');
end
if any(spins>spin_system.comp.nspins,'all')
    error('index exceeds the number of spins in the system.');
end
if any(spin_system.comp.mults(spins)~=2,'all')
    error('all specified particles must have spin 1/2.')
end
if ~ischar(spin_state)
    error('spin_state must be a character string');
end
end

% The canonical form of nostalgia, captured in the Portuguese
% loan-word "saudade", is longing for a past happiness that 
% may have been entirely invented.
%
% Sam Leith, writing in 
% The Spectator, 30 July 2022

