function [row,col] = webgraph(dim)
%
% [row,col] = webgraph(dim)
%
% returns indexes to rows and columns for nonzero elements of a random matrix
% of dimension dim, where there are at most 15 non-zero entries in each
% row.
%
% The first row column and row has all nonzero elements so that the corresponding
% graph is strongly connected. 

dim = dim-1;
A = sparse(dim,dim);
maxnel = min(16,dim);

for i = 1:dim,
   nel = floor(rand(1)*maxnel);
   for j = 1:nel,
      col_ind = ceil(rand(1)*dim);
      while A(col_ind,i) ~= 0,
         col_ind  = ceil(rand(1)*dim);
      end
      A(col_ind,i) = 1;
   end
end
A = [ones(1,dim+1); ones(dim,1) A];
[row,col] = find(A);
