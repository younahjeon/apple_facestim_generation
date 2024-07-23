function Lighting(varargin)

% Adjust lighting and size of the mesh

% Lighting(cameraViewAngle,cameraTarget,cameraPos,lightingAngle,lightingIntensity)

if isempty(varargin)
    cameraTarget = [0,0,0];
    cameraPos = [0,0,0];
    cameraViewAngle = 7.5;
    viewAngle = 0; 
    lightingAngle = {[45-rad2deg(viewAngle),25],[-45-rad2deg(viewAngle),25]};
    lightingIntensity = {[0.7,0.7,0.7],[0.7,0.7,0.7]};
else 
    cameraViewAngle = varargin{1};
    cameraTarget = varargin{2};
    cameraPos = varargin{3};
    lightingAngle = varargin{4};
    lightingAngle = cellfun(@rad2deg,lightingAngle,'UniformOutput',false); 
    lightingIntensity = varargin{5};
end

set(gca,'CameraViewAngle',cameraViewAngle,'CameraTarget',cameraTarget,'CameraPosition',cameraPos);
axis equal
view (0,90)

for i = 1:size(lightingAngle,2)
    h = camlight(lightingAngle{i}(1),lightingAngle{i}(2));
    lighting gouraud;
    set(h,'color',lightingIntensity{i});
end
material([.3 .8 .1 10 1]);
set(gcf,'color','w'); 
axis off

end
