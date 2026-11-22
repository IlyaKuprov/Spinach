% Independent row-mask reference for magnetisation flux. Syntax:
%
%                  F=flux_row_mask(s,flux_rate)
%
% Parameters:
%
%    s         - based small spherical-tensor spin system with a flux type
%
%    flux_rate - N-by-N input rate matrix, in Hz; source rows and destination
%                columns are ordered by global spin index
%
% Outputs:
%
%    F         - the reference flux generator in the direct-sum basis
%
% ilya.kuprov@weizmann.ac.il
% talos@spindynamics.org

function F=flux_row_mask(s,flux_rate)

% Check consistency
grumble(s,flux_rate);

% Embed local descriptors only in this small independent regression reference
b=zeros(s.bas.offsets(end),s.comp.nspins);
for n=1:numel(s.chem.parts)
    rows=(s.bas.offsets(n)+1):s.bas.offsets(n+1);
    b(rows,s.chem.parts{n})=s.bas.basis{n};
end
orders=sum(logical(b),2); F=sparse(size(b,1),size(b,1));
[from,to,rates]=find(flux_rate);
for n=1:numel(rates)
    for order=1:2
        mask=(orders==1);
        if order==2, mask=(orders>1); end
        source=mask&(b(:,from(n))~=0);
        destin=mask&(b(:,to(n))~=0);
        if (order==1)||strcmp(s.rlx.flux_type,'intramolecular')
            stationary=source&destin;
            source=xor(source,stationary); destin=xor(destin,stationary);
            if ~isequal(b(source,from(n)),b(destin,to(n)))
                error('reference tensor labels do not match.');
            end
            F=F+rates(n)*sparse(find(destin),find(source),ones(nnz(source),1),size(F,1),size(F,2));
        end
        F=F-rates(n)*spdiags(double(source),0,size(F,1),size(F,2));
    end
end

end

% Consistency enforcement
function grumble(s,flux_rate)
if (~isstruct(s))||(~isfield(s,'bas'))
    error('s must be a Spinach spin system with a basis.');
end
if (~isnumeric(flux_rate))||(~isreal(flux_rate))||...
   (~isequal(size(flux_rate),[s.comp.nspins s.comp.nspins]))
    error('flux_rate must be a real square rate matrix with one row per spin.');
end
end


