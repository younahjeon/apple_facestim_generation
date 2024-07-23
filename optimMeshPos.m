function[p_new, optim_param] = optimMeshPos(p,param_init,lb, ub)
% part vertices
% initial parameters
% parameter lower bound
% parameter upper bound
% refind of choice. different for parts 
% condition specifies which optimization cost function to use (see
% part2baseDist.m)

optim_param= fmincon(@part2baseDist,param_init,[],[],[],[],lb,ub);
T = eye(4);
T(:,4) = [optim_param(1);optim_param(2);optim_param(3);1];
p = [p,ones(size(p,1),1)];

pfinal = T*p';

p_new = pfinal(1:3,:)';
end


