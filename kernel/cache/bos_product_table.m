% Structure coefficient tables for the associative envelopes of trun-
% cated Weyl algebras. Bosonic monomials in boson_mono(nlevels) are
% orthognalised by boson_mono_ortho(nlevels)and then used to calculate 
% the structure coeffcients. 
% 
% Syntax:
%
% [product_table_left,product_table_right]=bos_product_table(nlevels)
% 
% Parameters:
%
%          nlevels - number of bosonic ladder population levels 
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
% Note: bosonic monomials have been orthogonalised; the overlap matrix 
%       is the identity matrix upto machine precision.
%
% sarbojoy.das@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bos_product_table.m>

function [product_table_left,product_table_right]= bos_product_table(nlevels)

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
    
    % Get orthogonal bosonic monomials
    B = boson_mono_ortho(nlevels);

    % Precompute norms
    norms=zeros(nlevels^2,1);
    for n=1:nlevels^2
        norms(n)=sqrt(hdot(B{n},B{n}));
    end
    
    % Preallocate IST product tables
    product_table_left= zeros(nlevels^2,nlevels^2,nlevels^2);
    product_table_right=zeros(nlevels^2,nlevels^2,nlevels^2);

    % Populate product tables
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

    try % Try to save a cache record, but do not insist
        save(table_file,'product_table_left',...
                        'product_table_right','-v7.3'); drawnow;
    catch
        warning('Spinach directory appears to be write-protected');
    end

end

end

% Consistency enforcement
function grumble(nlevels)
if (~isnumeric(nlevels))||(~isreal(nlevels))||...
   (~isscalar(nlevels))||(nlevels<1)||(mod(nlevels,1)~=0)
    error('nlevels must be a positive integer.');
end
end

% Rage Against the Machine [a rock band] never specified
% what type of machine they were furious with, but I rec-
% kon it was probably a printer.
%
% John Moynes