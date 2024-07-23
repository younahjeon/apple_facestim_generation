function [allrows,idx] = findRows(mat,val,varargin)
% find rows of mat that contain val regardless of the column position of
% val
% varargin
% 'and' : returns rows of mat that contain all of the elements in val
% 'or' : returns any rows of mat that have elements in val. these rows can
% have other numbers that don't exist in val 
% 'or+' : returns any rows of mat that have a subset of val 
% when length(val) = 1, the function doesn't need string input because
% 'and','or','or+' will return the same rows

% default varargin is 'and'

% ex) mat = [1,2,3;
%             4,7,8;
%              7,8,9];

%     val = [4,7,8,10];

%     findRows(mat,val,'and') returns []
%     findRows(mat,val,'or') returns second and third row of mat
%     findRows(mat,val,'or+') returns second row of mat

szmat = size(mat,1);

if szmat ~=0
idx = false(szmat,1);
andor = varargin; 

for i = 1:size(mat,1)
    if strcmp(andor,'and') 
        idx(i) = all(ismember(val,mat(i,:)));
    elseif strcmp(andor,'or')
        idx(i) = any(ismember(val,mat(i,:)));
    elseif strcmp(andor, 'or+') 
        idx(i) = all(ismember(mat(i,:),val));                        
    else % length(val) = 1
        idx(i) = all(ismember(val,mat(i,:)));
    end
end

allrows = mat(idx,:); 
else % if mat is empty 
    fprintf('your matrix should have size greater than 0') 
end

end
