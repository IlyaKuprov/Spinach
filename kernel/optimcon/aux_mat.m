% Builds auxiliary matrices for the calculation of the directional 
% derivatives of the trapezium product quadrature propagator:
%
%           expm(-1i*((HL+HR)/2+(1i*dt/12)*[HL,HR])*dt)
% 
% with respect to control coefficients in the evolution generators
% HL and HR on the left and the right edge of the interval. The de-
% rivatives are calculated using Eq 16 of Goodwin and Kuprov:
%
%                 https://doi.org/10.1063/1.4928978
%
% Syntax:
%
%   [auxm_l,auxm_r]=aux_mat(drifts,controls,cc_comm_idx,...
%                                  cc_comm,dt,cL,cR,k,j)
%
% Parameters:
%
%     drifts - a cell array of two matrices containing drift 
%              generators at the left (first element) and the
%              right (second element) edge of the interval
%
%   controls - a cell array of K control generators
%
%cc_comm_idx - a KxK matrix of logicals indicating non-zero 
%              commutation of controls
%
%    cc_comm - a KxK cell array control commutation relarions
%         
%         dt - interval duration, seconds
%
%         cL - control generator coefficients at the left 
%              edge of the interval
%
%         cR - control generator coefficients at the right
%              edge of the interval
%
%          k - the index of the generator inside controls 
%              array that the differentiation refers to
%
%          j - (optional) the index of the 2nd generator 
%              inside controls array that the differentiation 
%              refers to. Required for 3x3 block auxiliary 
%              matrices
%
% Outputs:
%
%     auxm_l - auxilary matrix for the derivative of the
%              interval propagator with respect to cL
%
%     auxm_r - auxilary matrix for the derivative of the
%              interval propagator with respect to cR
%
% david.goodwin@inano.au.dk
% u.rasulov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=aux_mat.m>

function [auxm_l,auxm_r]=aux_mat(drifts,controls,cc_comm_idx,...
                                        cc_comm,dt,cL,cR,k,j)

% If no second index, define as 0
if nargin<9, j=0; end

% Check consistency
grumble(drifts,controls,cc_comm_idx,cc_comm,dt,cL,cR,k,j);

% Build left and right generators
G{1}=drifts{1}; G{2}=drifts{2};
for n=1:numel(controls)
    G{1}=G{1}+cL(n)*controls{n};
    G{2}=G{2}+cR(n)*controls{n};
end

% Build the overall generator and a blank matrix
G=(G{1}+G{2})/2+1i*dt*(1/12)*(G{1}*G{2}-G{2}*G{1});
B=spalloc(size(G,1),size(G,2),0);

% Start right directional derivative of the generator
Hk_dir_R=(1/2)*controls{k}+1i*dt*(1/12)*(drifts{1}*controls{k}-...
                                         controls{k}*drifts{1});

% Start left directional derivative of the generator
Hk_dir_L=(1/2)*controls{k}+1i*dt*(1/12)*(controls{k}*drifts{2}-...
                                         drifts{2}*controls{k});

% Contributions from other controls
for n=1:numel(controls)

    % If the commutator is non-zero
    if ~cc_comm_idx(n,k)
        
        % Add right side contribution
        Hk_dir_R=Hk_dir_R+1i*dt*(1/12)*cL(n)*cc_comm{n,k};
    
        % Add left side contribution
        Hk_dir_L=Hk_dir_L+1i*dt*(1/12)*cR(n)*cc_comm{k,n};

    end

end

% For 3x3 block matrices
if j~=0
    
    % Start right directional derivative of the generator
    Hj_dir_R=(1/2)*controls{j}+1i*dt*(1/12)*(drifts{1}*controls{j}-...
                                             controls{j}*drifts{1});
    
    % Start left directional derivative of the generator
    Hj_dir_L=(1/2)*controls{j}+1i*dt*(1/12)*(controls{j}*drifts{2}-...
                                             drifts{2}*controls{j});
    
    % Contributions from other controls
    for n=1:numel(controls)

        % If the commutator is non-zero
        if ~cc_comm_idx(n,j)

            % Add right side contribution
            Hj_dir_R=Hj_dir_R+1i*dt*(1/12)*cL(n)*cc_comm{n,j};

            % Add left side contribution
            Hj_dir_L=Hj_dir_L+1i*dt*(1/12)*cR(n)*cc_comm{j,n};

        end

    end
    
end

% Build auxiliary matrices
if j==0

    % 2x2 block matrix
    auxm_l=[G, Hk_dir_L; B, G];
    auxm_r=[G, Hk_dir_R; B, G];

else
    
    % 3x3 block matrix
    auxm_l=[G, Hk_dir_L, B; B, G, Hj_dir_L; B, B, G];
    auxm_r=[G, Hk_dir_R, B; B, G, Hj_dir_R; B, B, G];

end

end

% Consistency enforcement
function grumble(drifts,controls,cc_comm_idx,cc_comm,dt,cL,cR,k,j)
if ~iscell(drifts)
    error('drifts must be a cell array of matrices.');
end
for n=1:numel(drifts)     
    if (~isnumeric(drifts{n}))||(size(drifts{n},1)~=size(drifts{n},2))
        error('all elements of drifts cell array must be square matrices.');
    end
end
if ~iscell(controls)
    error('controls must be a cell array of square matrices.');
end
for n=1:numel(controls)
    if (~isnumeric(controls{n}))||...
       (size(controls{n},1)~=size(controls{n},2))||...
       (size(controls{n},1)~=size(drifts{1},1))
        error('control generators must have the same size as drift generators.');
    end
end
if (~islogical(cc_comm_idx))||...
   (size(cc_comm_idx,1)~=size(cc_comm_idx,2))||...
   (size(cc_comm_idx,1)~=numel(controls))
    error('cc_comm_idx must be a square matrix of logicals, with length equal to the number of controls.');
end
if ~iscell(cc_comm)
    error('cc_comm must be a cell array of square matrices.');
end
for n=1:numel(cc_comm)
    if (~isnumeric(cc_comm{n}))||...
       (size(cc_comm{n},1)~=size(cc_comm{n},2))||...
       (size(cc_comm{n},1)~=size(drifts{1},1))
        error('cc_comm generators must have the same size as drift generators.');
    end
end
if (~isnumeric(dt))||(~isreal(dt))||(~isscalar(dt))||(dt<=0)
    error('dt must be a positive real scalar.');
end
if (~isnumeric(cL))||(~isreal(cL))||(~isnumeric(cR))||(~isreal(cR))
    error('cL and cR must be arrays of real numbers.');
end
if (~isnumeric(k))||(~isreal(k))||(~isscalar(k))||(k<=0)||(mod(k,1)~=0)
    error('k must be a positive real integer.');
end
if (~isnumeric(j))||(~isreal(j))||(~isscalar(j))||(j<0)||(mod(j,1)~=0)
    error('j must be a non-negative real integer.');
end
end

% Rutherford to Eddington (who wondered at table whether we should ever
% come to know electrons as more than mental concepts): "Not exist? Not
% exist?! Why, I can see the little buggers as plain as I can see that 
% spoon in front of me!"

