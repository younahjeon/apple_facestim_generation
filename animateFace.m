function animateFace(faceData,varargin)
% Animate all sessions
% 
load('TriangleIndices.mat');
figure;
if nargin>1
    if ~strcmp(varargin{1},'video')
        error('wrong input')
    end
    v = VideoWriter('movie.avi');
    open(v)
end

for ns = 1:length(faceData.vertices)
    ns
    vertices = faceData.vertices{ns};
   
    trimesh(TriangleIndices',vertices(:,1),vertices(:,2),vertices(:,3),...
         'FaceColor',[0.65,0.65,0.65],'EdgeColor','none');
    
    cameraViewAngle = 10; % field of view (0 to 180) The greater the number, smaller objects will show more
    cameraTarget = [32 10.785 14.3812]; % how far your target is from the camera. determines the degree of visibility
    cameraPos = [0 0 0];
    set(gca,'CameraViewAngle',cameraViewAngle,'CameraTarget',cameraTarget,'CameraPosition',cameraPos);
    axis equal
    view (0,90)
    
    h1 = camlight(45,25); lighting gouraud;
    set(h1,'color',[.7 .7 .7]);
    h3 = camlight('headlight'); lighting gouraud
    set(h3,'color',[0.2 0.2 0.2]);
    h2=camlight(-45,25); lighting gouraud
    set(h2,'color',[0.7 0.7 0.7]);
    material([.3 .8 .1 10 1]);
    set(gcf,'color','w');
    set(gcf,'menubar','none');
    axis off
    if i ==1
    set(gcf,'position',[1 1 1.15 1.15].*get(gcf,'position'));
    view (0,90)
    end
    title(sprintf('%s frame = %d',faceData.ID,ns));
    drawnow 
    pause(0.01)
    
    if length(varargin)>1
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end
end
