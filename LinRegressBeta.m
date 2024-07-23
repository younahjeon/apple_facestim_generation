function [beta] = LinRegressBeta(Y,X,ridgeweight)
%% Finding Linear Mapping between Mesh and BlendShapes


% Y = Xb
% Y is mesh vertices.
% size = number of Frames x number of Vertices (nf x 3660)
%
% X is blendshape coeffiecients for a specific mesh.
% size = number of Frames x number of BlendShapes (nf x 53)
% a dummy variable is needed, so X should have 1's at 53rd column.
%
% b is betas(weights), or how much each blendshape coefficient influences a
% specific vertex. Each row represents how much a single blendshape
% coefficient distorts a mesh
% size = number of BlendShapes x number of Verties (53 x 3660)

% ridgeweight determines the strength of L2 norm. Ridge regression solves
% the multicollinearity problem which arises when independent variables are
% actually near-linearly related. With enough samples (number of frames), you can
% avoid this problem. 
    beta = pinv(X'*X+ ridgeweight*eye(size(X,2)))*X'*Y;
end




