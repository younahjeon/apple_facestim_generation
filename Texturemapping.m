% map texture onto the mesh
fname = 'faceData_example.txt';
faceData = faceData_readLog(fname); 

nVertices = 1220;

frame_to_use = 1; 

%transform the vertex to world coordinate
v = faceData.vertices{frame_to_use};
nv = [v, ones(nVertices,1)];
ct =faceData.ct{frame_to_use};

ft = faceData.ft{frame_to_use};
cV = ft * nv';

% from world to camera
cV = pinv(ct) * cV;
cV = cV(1:3,:)';

% camera to image 
pt= faceData.pt{frame_to_use};
nv = [cV, ones(nVertices,1)];
projectedCoord=  pt * nv'; 
projectedCoord= projectedCoord(1:3,:)./projectedCoord(4,:); % normalize
projectedCoord = projectedCoord'; 

% change to uv coordinate
x = projectedCoord(:,1) *0.5 + 0.5;
y = projectedCoord(:,2) * 0.5 + 0.5; 
y = 1 - y; % image coordinates are y-flipped 

% to pixels 
pixelCoord = zeros(nVertices,1);
pixelCoord(:,1) = x*1440;
pixelCoord(:,2) = y *1080;

Options = struct;
I= imread('0000.jpg');
load('TriangleIndices.mat')
patcht(TriangleIndices', v, TriangleIndices',flip(pixelCoord,2),I, Options)
