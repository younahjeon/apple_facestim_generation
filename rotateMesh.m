function p_new = rotateMesh(p, rotType,varargin)
% rotates a set of 3d points 
% p is a 3-dimensional array 
% rotType can be axis-angle, for which the user provides the rotation axis
% and angle
% rotType can be vertexnormal, for which the user provides the starting and
% ending vertex normal. 
global partIndices

    if strcmp(rotType,'axis-angle')
        rotAxis = varargin{1};
        rotAngle = varargin{2};
    elseif strcmp(rotType,'vertexnormal')
        refVN = varargin{1};
        targetVN = varargin{2};
        rotAxis = cross(refVN,targetVN); 
        normalVN = cross(refVN,targetVN); % normal of the refVN and targetVN
        rotAngle = acos(dot(refVN,targetVN)); 
        if dot(rotAxis,normalVN) <0  
            rotAngle = -rotAngle;
        end
    end
    
    R = vrrotvec2mat([rotAxis,rotAngle]);
    R = [[R;0,0,0],[0;0;0;1]];
    
    center_delta = mean(p);
    p = p-center_delta;
    p = [p,ones(size(p,1),1)];
    pfinal = R*p';

    p_new = pfinal(1:3,:)';
    
    figure; 
    h1 = trimesh(partIndices,p(:,1),p(:,2),p(:,3),'FaceVertexCData',[1,1,1],'EdgeColor',[0,0,0],'FaceColor','flat');
    hold on
    h2 = trimesh(partIndices,p_new(:,1),p_new(:,2),p_new(:,3),'FaceVertexCData',[1,1,1],'EdgeColor',[0.5,0.5,0.5],'FaceColor','flat');
    
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend([h1,h2],{'original','updated'})
    xlim([min(p(:,1))-0.05, max(p(:,1))+0.05])
    ylim([min(p(:,2))-0.05, max(p(:,2))+0.05])
    zlim([min(p(:,3))-0.05, max(p(:,3))+0.05])


end