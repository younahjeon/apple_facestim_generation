function [targetCoord targetVN] = targetSelection(F,V)  

global baseTR baseVN

h = figure;
trimesh(F,V(:,1),V(:,2),V(:,3),...
				'FaceVertexCData',[1,1,1],'EdgeColor',[0,0,0],'FaceColor','flat',...
            'FaceAlpha',1); 
        cameraViewAngle = 10; % field of view (0 to 180) The greater the number, smaller objects will show more
cameraTarget = [32 10.785 14.3812]; % how far your target is from the camera. determines the degree of visibility
cameraPos = [0 0 0];
set(gca,'CameraViewAngle',cameraViewAngle,'CameraTarget',cameraTarget,'CameraPosition',cameraPos);
axis equal
view (0,90) 

dcm_obj = datacursormode(h); 
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on')

disp('Click a point to display a data tip, then press Return.')

pause 
info_struct = getCursorInfo(dcm_obj);
targetCoord = info_struct(1).Position;
centerOccluded = nearestNeighbor(baseTR,targetCoord);
targetVN = baseVN(centerOccluded,:); 

close(h)
end