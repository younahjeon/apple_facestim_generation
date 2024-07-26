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

%% add expression to a textured mesh

newv = changeMeshbyEmotion('smile',v, TriangleIndices');

figure;
subplot(1,2,1);
patcht(TriangleIndices', v, TriangleIndices',flip(pixelCoord,2),I, Options)
cameraViewAngle = 10; % field of view (0 to 180) The greater the number, smaller objects will show more
cameraTarget = [32 10.785 14.3812]; % how far your target is from the camera. determines the degree of visibility
cameraPos = [0 0 0];
set(gca,'CameraViewAngle',cameraViewAngle,'CameraTarget',cameraTarget,'CameraPosition',cameraPos);
axis equal
view (0,90)

%h1 = camlight(45,25); lighting gouraud;
%set(h1,'color',[.7 .7 .7]);
%h3 = camlight('headlight'); lighting gouraud
%set(h3,'color',[0.2 0.2 0.2]);
%h2=camlight(-45,25); lighting gouraud
%set(h2,'color',[0.7 0.7 0.7]);
%material([.3 .8 .1 10 1]);
set(gcf,'color','w');
set(gcf,'menubar','none');
axis off
title('original')

subplot(1,2,2);
patcht(TriangleIndices', newv, TriangleIndices',flip(pixelCoord,2),I, Options)
cameraViewAngle = 10; % field of view (0 to 180) The greater the number, smaller objects will show more
cameraTarget = [32 10.785 14.3812]; % how far your target is from the camera. determines the degree of visibility
cameraPos = [0 0 0];
set(gca,'CameraViewAngle',cameraViewAngle,'CameraTarget',cameraTarget,'CameraPosition',cameraPos);
axis equal
view (0,90)

%h1 = camlight(45,25); lighting gouraud;
%set(h1,'color',[.7 .7 .7]);
%h3 = camlight('headlight'); lighting gouraud
%set(h3,'color',[0.2 0.2 0.2]);
%h2=camlight(-45,25); lighting gouraud
%set(h2,'color',[0.7 0.7 0.7]);
%material([.3 .8 .1 10 1]);
set(gcf,'color','w');
set(gcf,'menubar','none');
axis off
title('original + smile Blendshapes')

savefig('textured_mesh_expresison.png')
