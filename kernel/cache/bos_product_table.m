% Structure coefficient tables for the associative envelopes of trun-
% cated Weyl algebras. Syntax:
%
% [product_table_left,product_table_right]=bos_product_table(nlevels)
%
% The input parameter is the number of energy levels in the truncated 
% bosonic mode. Output contains the structure coefficients in the fol-
% lowing conventions (B{m} and B{k} are normalised):
%
%            B{n}*B{m}=...+product_table_left(n,m,k)*B{k}+...
%
%            B{m}*B{n}=...+product_table_right(n,m,k)*B{k}+...
%
% corresponding to the expansion of the left and the right multiplica-
% tive action by B{n} on B{m} as given in Eq 7.18 of the first edition
% of IK's book (normalisation is missing in the book, that's a typo).
% Numbering translation between single and double index is given by 
% kq2lin and lin2kq functions.
%
% Note: these are expensive tables, disk cache is used automatically.
%
% Note: unlike irreducible spherical tensors, bosonic monomials are
%       not orthogonal; the overlap matrix is built and used below.
%
% sarbojoy.das@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bos_product_table.m>

function [product_table_left,product_table_right]=bos_product_table(nlevels)

% Check consistency
grumble(nlevels);

% Generate cache record name
own_path=mfilename('fullpath');
own_path=own_path(1:(end-17));
table_file=[own_path 'bos_product_table_' ...
            num2str(nlevels) '.mat'];

% Check the cache
if exist(table_file,'file')
    
    % Lift data from the cache if the file is already available
    load(table_file,'product_table_left','product_table_right');
    
else
    
    % Get bosonic monomials
    B=boson_mono(nlevels);

    % Precompute norms
    norms=zeros(nlevels^2,1);
    for n=1:nlevels^2
        norms(n)=sqrt(hdot(B{n},B{n}));
    end

    % Get the overlap matrix
    S=zeros(nlevels^2,nlevels^2);
    for m=1:nlevels^2
        for k=1:nlevels^2

            % Carefully tiptoe around extreme norms
            S(m,k)=hdot(B{m}/norms(m),B{k}/norms(k));

        end
    end
    
    % Preallocate the arrays
    product_table_left= zeros(nlevels^2,nlevels^2,nlevels^2);
    product_table_right=zeros(nlevels^2,nlevels^2,nlevels^2);
    
    % Get the structure coefficients
    for k=1:nlevels^2
        for m=1:nlevels^2
            for n=1:nlevels^2

                % Left product action: carefully tiptoe around extreme norms
                product_table_left(n,m,k)= norms(n)*hdot((B{k}/norms(k)),...
                                                         (B{n}/norms(n))*...
                                                         (B{m}/norms(m)));

                % Right product action: carefully tiptoe around extreme norms
                product_table_right(n,m,k)=norms(n)*hdot((B{k}/norms(k)),...
                                                         (B{m}/norms(m))*...
                                                         (B{n}/norms(n)));

            end
        end
    end

    % Apply the overlap matrix
    product_table_left= reshape(product_table_left, [nlevels^4 nlevels^2]);
    product_table_right=reshape(product_table_right,[nlevels^4 nlevels^2]);
    product_table_left= product_table_left /S';
    product_table_right=product_table_right/S';
    product_table_left= reshape(product_table_left, [nlevels^2 nlevels^2 nlevels^2]);
    product_table_right=reshape(product_table_right,[nlevels^2 nlevels^2 nlevels^2]);

    try % Try to save a cache record, but don't insist
        save(table_file,'product_table_left',...
                        'product_table_right'); drawnow;
    catch
        warning('Spinach installation appears to be write-protected');
    end

end

end

% Consistency enforcement
function grumble(nlevels)
if (~isnumeric(nlevels))||(~isreal(nlevels))||...
   (~isscalar(nlevels))||(nlevels<3)||(mod(nlevels,1)~=0)
    error('mult must be a real integer greater than 2.');
end
end

% Rage Against the Machine [a rock band] never specified
% what type of machine they were furious with, but I rec-
% kon it was probably a printer.
%
% John Moynes

