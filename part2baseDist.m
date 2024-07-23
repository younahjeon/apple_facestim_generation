function sumDist = part2baseDist(params)
% calculate the mean squared error of the distance between the base face
% mesh and the part mesh
% for now, we only allow translations for optimization.(no rotations)

global condition refInd_optim partB targetCoord baseTR partV

transx = params(1);
transy = params(2);
transz = params(3);
%   thetax = params(4);
%   thetay = params(5);
%   thetaz = params(6);
%
%   Rx = [1,0,0 0;...
%     0,cos(thetax),-1*sin(thetax),0;...
%     0,sin(thetax),cos(thetax),0;...
%     0,0,0,1];
%
%  Ry = [cos(thetay),0,sin(thetay),0;...
%     0,1,0,0;...
%     -sin(thetay),0,cos(thetay),0;...
%     0,0,0,1];
%
%  Rz = [cos(thetaz),-sin(thetaz),0,0;...
%     sin(thetaz),cos(thetaz),0,0;...
%     0,0,1,0;
%     0,0,0,1];

T = eye(4);
T(:,4) = [transx;transy;transz;1];
% Tfinal = T*Rx*Ry*Rz;
Tfinal = T;
p= [partV,ones(size(partV,1),1)];
pfinal = Tfinal*p';

if strcmp(condition,'optim1')
    minDistArr = pdist2(pfinal(1:3,refInd_optim)',targetCoord);
    sumDist =sum(minDistArr);
elseif strcmp(condition,'optim2')
    % optimizes for the distance between part boundary vertices and the base 
    p_boundary = pfinal(1:3,partB);
    [~,d] = nearestNeighbor(baseTR,p_boundary');
    minDistArr = d;
    sumDist = sum(minDistArr);
else
    error('wrong condition input')
end

end