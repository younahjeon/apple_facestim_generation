function animateFace_with_pause(faceData,varargin)
% Animate all sessions
global TriangleIndices;
load('TriangleIndices.mat');
figure;

% Initialize pause flag
paused = false;

% Add a button to pause/resume animation
uicontrol('Style', 'pushbutton', 'String', 'Pause/Resume',...
    'Position', [20 20 100 30],...
    'Callback', @pause_resume);

if nargin>1
    if ~strcmp(varargin{1},'video')
        error('wrong input')
    end
    v = VideoWriter('movie.avi');
    open(v)
end

for ns = 1:length(faceData.vertices)
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
    if ns ==1
        set(gcf,'position',[1 1 1.15 1.15].*get(gcf,'position'));
        view (0,90)
    end
    title(sprintf('frame = %d',ns));
    drawnow
    
    % Check if paused and wait until unpaused
    while paused
        pause(0.1);
    end
    
    if nargin>1 && strcmp(varargin{1},'video')
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    
    pause(0.01)
end

if nargin>1 && strcmp(varargin{1},'video')
    close(v)
end

    % Function to toggle pause/resume
    function pause_resume(~, ~)
        paused = ~paused; % Toggle pause state
    end
end