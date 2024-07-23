function [new_index] = reIndexing(old_index,varargin)
% reindexes 
% default is reindexing starts from 1 and the max index is the number
% of unique elements in old_index.
% if there are two inputs, the second input should be the reference index
% from which new_index is created. 
% in this case, length(reference index) should be equal to the number of
% unique elements in old_index 


lin_ind = old_index(:); 


if nargin==1
    unique_ind = unique(old_index);
else
    unique_ind = varargin{1};
end

new_index = zeros(length(lin_ind),1); 
for i = 1:length(lin_ind)
    if nargin ==1 
    new_index(i) = find((lin_ind(i)-unique_ind)==0); 
    else
         new_index(i) = unique_ind(lin_ind(i));
    end
end
new_index = reshape(new_index,size(old_index));
end


