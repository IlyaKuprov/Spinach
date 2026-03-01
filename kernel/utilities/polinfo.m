% Draws an ASCII diagram of a polyadic object. Syntax:
%
%                       polinfo(p)
%
% Parameters:
%
%       p   -   polyadic object
%
% Outputs:
%
%       an ASCII diagram to the console
%
% Note: polyadic objects can be huge, the code below 
%       avoids making memory copies.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=polyadic/polinfo.m>

function polinfo(p,level,label)

% Default settings
if ~exist('level','var'), level=0; end
if ~exist('label','var'), label='polyadic'; end

% Check consistency
grumble(p,level,label);

% Set indentation level
indent=repmat(' ',1,4*level);

% Get the dimensions
[nrows,ncols]=size(p);

% Print label and size
fprintf('%s%s [%dx%d]\n',indent,label,nrows,ncols);

% Print prefixes
if isa(p,'polyadic')
    if ~isempty(p.prefix)
        fprintf('%s    prefix:\n',indent);
    end
    for n=1:numel(p.prefix)
        if isa(p.prefix{n},'polyadic')
            polinfo(p.prefix{n},level+2,sprintf('polyad %d',n));
        elseif isa(p.prefix{n},'opium')
            [nrows,ncols]=size(p.prefix{n});
            fprintf('%s        opium  %d [%dx%d]\n',indent,n,nrows,ncols);
        else
            [nrows,ncols]=size(p.prefix{n});
            fprintf('%s        matrix %d [%dx%d]\n',indent,n,nrows,ncols);
        end
    end
end

% Print kron terms
for n=1:numel(p.cores)
    fprintf('%s    kron %d\n',indent,n);
    for k=1:numel(p.cores{n})
        if isa(p.cores{n}{k},'polyadic')
            polinfo(p.cores{n}{k},level+2,sprintf('polyad %d',k));
        elseif isa(p.cores{n}{k},'opium')
            [nrows,ncols]=size(p.cores{n}{k});
            fprintf('%s        opium  %d [%dx%d]\n',indent,k,nrows,ncols);
        else
            [nrows,ncols]=size(p.cores{n}{k});
            fprintf('%s        matrix %d [%dx%d]\n',indent,k,nrows,ncols);
        end
    end
end

% Print suffixes
if isa(p,'polyadic')
    if ~isempty(p.suffix)
        fprintf('%s    suffix:\n',indent);
    end
    for n=1:numel(p.suffix)
        if isa(p.suffix{n},'polyadic')
            polinfo(p.suffix{n},level+2,sprintf('polyad %d',n));
        elseif isa(p.suffix{n},'opium')
            [nrows,ncols]=size(p.suffix{n});
            fprintf('%s    opium  %d [%dx%d]\n',indent,n,nrows,ncols);
        else
            [nrows,ncols]=size(p.suffix{n});
            fprintf('%s    matrix %d [%dx%d]\n',indent,n,nrows,ncols);
        end
    end
end

end

% Consistency enforcement
function grumble(p,level,label)
if ~isa(p,'polyadic')
    error('p must be a polyadic object.');
end
if (~isnumeric(level))||(~isreal(level))||(~isscalar(level))||...
   (level<0)||(mod(level,1)~=0)
    error('level must be a non-negative integer.');
end
if ~ischar(label)
    error('label must be a character string.');
end
end

% Ignoring bad ideas doesn't make them go away; they will 
% still eat up funding. [...] Killing ideas is a necessary
% part of science. Think of it as a community service.
%
% Sabine Hossenfelder, "Lost in Math"

