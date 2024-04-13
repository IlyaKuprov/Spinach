% Sums buffered tensor trains in a single tensor train using
% AMEn algorithm. Syntax:
%
%                      y=amensum(x,tol,opts)
%
% Parameters: 
%
%    x   - ttclass with buffered rank-one tensors
%
%    tol - relative tolerance parameter, e.g. 1e-10
%
% The opts field is optional:
%
%    opts.max_swp - maximum number of iterations
%
%    opts.init_guess_rank - rank of the initial guess
%
%    opts.enrichment_rank - rank of the enrichment
%
%    opts.verb            - verbosity switch
%
% Outputs: 
%
%    y   - ttclass with a single tensor train, such
%          that |x-y|<tol*|x| in Frobenius norm
%
% d.savosyanov@soton.ac.uk
% sergey.v.dolgov@gmail.com
%
% <https://spindynamics.org/wiki/index.php?title=ttclass/amensum.m>

function y=amensum(x,tol,opts)

% Process options
if ~exist('opts','var'), opts = struct; end
if ~isfield(opts, 'max_swp');         opts.max_swp=100;       end
if ~isfield(opts, 'init_guess_rank'); opts.init_guess_rank=2; end
if ~isfield(opts, 'enrichment_rank'); opts.enrichment_rank=4; end
if ~isfield(opts, 'verb');            opts.verb=0;            end

% Read tensor train sizes and ranks
[d,N]=size(x.cores); sz=x.sizes; rnk=x.ranks;

% Proceed a CP format and a sum of TT formats separately
if all(all(rnk==1))
    
    % Convert ttclass to CP format
    X=cell(d,1);
    for k=1:d
        X{k,1}=zeros(sz(k,1), sz(k,2),N);
        for i=1:N
            X{k,1}(:,:,i)=reshape(x.cores{k,i}, [sz(k,1),sz(k,2),1]);
        end
    end
    
    % Generate a random initial guess
    [y,~]=ttort(rand(x, opts.init_guess_rank),-1);
    
    % Generate a random enrichment
    [z,~]=ttort(rand(x, opts.enrichment_rank),-1);
    
    % Precompute interfaces of projections Y'X Z'X and Z'Y
    yx = cell(d+1,1); yx{1}=ones(1,N); yx{d+1}=ones(N,1);
    zx = cell(d+1,1); zx{1}=ones(1,N); zx{d+1}=ones(N,1);
    zy = cell(d+1,1); zy{1}=1; zy{d+1}=1;
    nrm = ones(d,1);
    for k=d:-1:2
        yx{k}=step_tt_by_cp(yx{k+1},y{k,1},X{k,1},-1);
        zx{k}=step_tt_by_cp(zx{k+1},z{k,1},X{k,1},-1);
        zy{k}=step_tt_by_tt(zy{k+1},z{k,1},y{k,1},-1);
        nrm(k)=norm(yx{k},'fro');
        if(nrm(k)>0)
            yx{k}=yx{k}/nrm(k);
            zx{k}=zx{k}/nrm(k);
        end
    end
    
    % Set initial conditions
    flipped=false; satisfied=false; iter=0;
    
    % Main iteration cycle
    while ~satisfied
        
        % Increment the counter
        iter=iter+1;
        
        % Read current ranks of approximation
        ry=y.ranks; rz=z.ranks; sz=x.sizes;
        
        % The difference between two consequitive iterations
        largest_stepsize=0;
        
        % Loop over the cores
        for k=1:d
            
            % Solve local approximation problem
            ynew=local_cp_to_tt(yx{k},yx{k+1},X{k,1},x.coeff);
            nrm(k)=norm(ynew(:),'fro');
            if (nrm(k)>0)
                ynew=ynew/nrm(k);
            else
                nrm(k)=1;
            end
            
            % Check the convergence
            relative_stepsize=norm(ynew(:)-y{k,1}(:),2)/norm(y{k,1}(:),2);
            largest_stepsize=max(largest_stepsize, relative_stepsize);
            
            % Run the truncation
            if k<d
                
                % Truncate the updated core
                Y = reshape(ynew,[ry(k)*sz(k,1)*sz(k,2),ry(k+1)]);
                [U,S,V] = svd(Y, 'econ'); S = diag(S);
                rnew = frob_chop(S,tol*norm(S,2)/sqrt(d));
                U = U(:,1:rnew); V = diag(S(1:rnew))*V(:,1:rnew)'; 
                ynew = U*V;
                
            end
            
            % Prepare error update and enrichment
            if opts.enrichment_rank>0
                
                % Compute new error core after truncation
                % Project the truncated solution to the interface of Z
                ynew = reshape(ynew, [ry(k)*sz(k,1)*sz(k,2)*ry(k+1), 1]);
                ynew_enrich = reshape(ynew, [ry(k)*sz(k,1)*sz(k,2), ry(k+1)]);
                ynew_enrich = ynew_enrich*zy{k+1};
                
                % Save crys. We use it in the enrichment of the solution
                yznew = reshape(ynew_enrich, [ry(k), sz(k,1)*sz(k,2)*rz(k+1)]);
                yznew = zy{k}*yznew;
                yznew = reshape(yznew, [rz(k)*sz(k,1)*sz(k,2)*rz(k+1), 1]);
                znew = local_cp_to_tt(zx{k},zx{k+1},X{k,1},x.coeff);
                znew = znew(:)/nrm(k) - yznew;
                znew = reshape(znew, [rz(k)*sz(k,1)*sz(k,2), rz(k+1)]);
                [znew,~]=qr(znew, 0); rz(k+1) = size(znew, 2);
                znew = reshape(znew, [rz(k),sz(k,1),sz(k,2),rz(k+1)]);
                
                % Run the enrichment
                if k<d
                    
                    % Compute the enrichment core
                    enrichment = local_cp_to_tt(yx{k},zx{k+1},X{k,1},x.coeff);
                    enrichment = enrichment(:)/nrm(k) - ynew_enrich(:);
                    enrichment = reshape(enrichment, [ry(k)*sz(k,1)*sz(k,2), rz(k+1)]);
                    
                    % Enrich and orthogonalize the solution
                    U_enriched = [U, enrichment]; [U,R]=qr(U_enriched, 0);
                    R = R(:, 1:rnew); V = R*V; rnew = size(U,2);
                    
                end
                
            end
            
            % Run the updates
            if k<d
                
                % Save the new block U
                y.cores{k,1} = reshape(U, [ry(k),sz(k,1),sz(k,2),rnew]);
                
                % Push the non-orthogonal factor forth the chain                
                next_core = reshape(y{k+1,1}, [ry(k+1), sz(k+1,1)*sz(k+1,2)*ry(k+2)]);
                next_core = V*next_core;
                y.cores{k+1,1} = reshape(next_core, [rnew,sz(k+1,1),sz(k+1,2),ry(k+2)]);
                
                % Update the rank
                ry(k+1)=rnew;
                
                % Update the interfaces
                yx{k+1}=step_tt_by_cp(yx{k},y{k,1},X{k,1},+1);
                nrm(k)=norm(yx{k+1},'fro');
                if (nrm(k)>0)
                    yx{k+1}=yx{k+1}/nrm(k);
                end
                if opts.enrichment_rank>0
                    zx{k+1}=step_tt_by_cp(zx{k},znew,X{k,1},+1);
                    zy{k+1}=step_tt_by_tt(zy{k},znew,y{k,1},+1);
                    if (nrm(k)>0)
                        zx{k+1}=zx{k+1}/nrm(k);
                    end
                end
                
            end
            
        end
        
        % Reverse the trains
        X = X(d:-1:1); y=revert(y); z=revert(z); x=revert(x);
        
        % Reverse the projections
        yx(2:d) = yx(d:-1:2);
        zx(2:d) = zx(d:-1:2);
        zy(2:d) = zy(d:-1:2);
        for k=2:d
            yx{k} = yx{k}.';
            zx{k} = zx{k}.';
            zy{k} = zy{k}.';
        end
        nrm=nrm(d:-1:1);
        
        % Mark that we are currently having a flipped TT
        flipped = ~flipped;
        
        % And exit only if the direction is correct
        satisfied = (largest_stepsize<tol) && (iter<=opts.max_swp);
        
        % Report if necessary
        if opts.verb>0
            str=sprintf('iter=%d, max_dx=%3.3e, max_r=%d\n',iter, largest_stepsize, max(y.ranks));
            report(spin_system,str);
        end
        
    end
    
    % Output in the ttclass format
    if flipped
        y=revert(y);
    end
    nrm=exp(sum(log(nrm))/d);
    for k=1:d
        y.cores{k,1}=nrm*y.cores{k,1};
    end
    y.tolerance=(nrm^d)*tol;
    
else
    
    % sum of tensor trains
    error('TT not ready yet');

end

end

% CP to TT conversion
function ttcore = local_cp_to_tt(left_interface,right_interface,cp,coeff)

% Read sizes
[r1,N1]=size(left_interface);
[N2,r2]=size(right_interface);
[sz1,sz2,Ncp]=size(cp);
[~,Ncf]=size(coeff);
if ~(N1==N2 && N2==Ncp && Ncp==Ncf)
    error('CP ranks mismatch.');
end
N=Ncf;

% Contract left intefrace with CP data
left_interface=repmat(left_interface,[sz1*sz2,1]); % [r1*sz1*sz2, N]
cp=reshape(cp,[1,sz1*sz2*N]);
cp=repmat(cp,[r1,1]);                              % [r1,sz1*sz2*N]
cp=reshape(cp,[r1*sz1*sz2, N]);
left_interface=left_interface.*cp;
left_interface=reshape(left_interface,[r1*sz1*sz2,N]);

% Contract right interface with coeff vector
coeff=reshape(coeff,[N,1]);
coeff=repmat(coeff,[1,r2]);
right_interface=right_interface.*coeff;            % [N,r2]

% Contract left and right interfaces
ttcore=left_interface*right_interface;
ttcore=reshape(ttcore,[r1,sz1,sz2,r2]);

end

% TT by CP
function next_interface = step_tt_by_cp(interface, tt, cp, direction)

% Read sizes
[r1,sz1,sz2,r2] = size(tt);
[nn1,nn2,N] = size(cp);
if ~(sz1==nn1 && sz2==nn2)
    error('sizes mismatch.');
end

% Decide direction
switch direction
    
    case +1
        
        % Move one core to the right, compute next left interface
        tt = reshape(tt, [r1, sz1*sz2*r2]);
        interface = reshape(interface, [r1,N]);
        next_interface = tt'*interface;
        next_interface = reshape(next_interface, [sz1*sz2*r2, N]);
        cp = reshape(cp, [sz1*sz2,N]);
        cp = repmat(cp, [r2,1]);
        next_interface = next_interface.*cp;
        next_interface = reshape(next_interface, [sz1*sz2, r2*N]);
        next_interface = sum(next_interface, 1);
        next_interface = reshape(next_interface, [r2, N]);
        
    case -1
        
        % Move one core to the left, compute next right interface
        tt = reshape(tt, [r1*sz1*sz2, r2]);
        interface = reshape(interface, [N,r2]);
        next_interface = interface*tt';
        next_interface = reshape(next_interface, [N, r1*sz1*sz2]);
        cp = reshape(cp, [sz1*sz2,N]);
        cp = cp.';
        cp = repmat(cp, [r1,1]);
        cp = reshape(cp, [N, r1*sz1*sz2]);
        next_interface = next_interface.*cp;
        next_interface = reshape(next_interface, [N*r1, sz1*sz2]);
        next_interface = sum(next_interface, 2);
        next_interface = reshape(next_interface, [N, r1]);
        
    otherwise
        
        % Complain and bomb out
        error('unrecognized direction.');
        
end

end

% TT by TT
function next_interface = step_tt_by_tt(interface, ttx, tty, direction)

% Read sizes
[rx1,szx1,szx2,rx2] = size(ttx);
[ry1,szy1,szy2,ry2] = size(tty);
if ~(szx1==szy1 && szx2==szy2)
    error('sizes mismatch.');
end
sz1=szx1; sz2=szx2;

% Decide direction
switch direction
    
    case +1 
        
        % Move one core to the right, compute next left interface
        ttx = reshape(ttx, [rx1, sz1*sz2*rx2]);
        tty = reshape(tty, [ry1*sz1*sz2, ry2]);
        interface = reshape(interface, [rx1, ry1]);
        next_interface = ttx'*interface;
        next_interface = reshape(next_interface, [sz1*sz2, rx2*ry1]);
        next_interface = next_interface.';
        next_interface = reshape(next_interface, [rx2, ry1*sz1*sz2]);
        next_interface = next_interface*tty;
        next_interface = reshape(next_interface, [rx2, ry2]);
        
    case -1 
        
        % Move one core to the left, compute next right interface
        ttx = reshape(ttx, [rx1, sz1*sz2*rx2]);
        tty = reshape(tty, [ry1*sz1*sz2, ry2]);
        interface = reshape(interface, [ry2, rx2]);
        next_interface = tty*interface;
        next_interface = reshape(next_interface, [ry1, sz1*sz2*rx2]);
        next_interface = next_interface*ttx';
        next_interface = reshape(next_interface, [ry1, rx1]);
        
    otherwise
        
        % Complain and bomb out
        error('unrecognized direction.');
        
end

end

% The doctrine of equality! [...] But there is no more venomous
% poison in existence: for it appears to be preached by justice
% itself, when it is actually the end of justice.
%
% Friedrich Nietzsche

