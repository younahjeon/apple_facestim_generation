function [newMesh,optim_param]= MeshSurgery(t,part_string,smoothingArea,newMesh)
global condition partB partIndices...
    baseF baseV baseTR baseVN baseVA targetCoord partV

cameraTarget = [0,0,0];
cameraPos = [0,0,0];
cameraViewAngle = 7.5;
%% (1) Optimize part location

% Bring the part towards the target location by minimizing the distance
% between the reference point and the target point

param_init = [0,0,0]; % x, y, z coordinate 
lb = [-inf,-inf,-inf];
ub = [inf,inf,inf];
condition = 'optim1';

[partV, optim_param] = optimMeshPos(partV,param_init,lb, ub);

if strcmp(part_string, 'NOSE')
    param_init = [0,0,0];
    var_y = 0.005;
    lb = [-inf,-var_y,-inf];
    ub = [inf,var_y,inf];
    condition = 'optim2'; 
    [partV,optim_param] = optimMeshPos(partV,param_init,lb,ub); 
end
%% (2) Find the faces in the base mesh that will be occluded by the part
[centerOccluded,occludedI] = intersectMesh();
%% (3) Remove occluded Faces
newPTvert = baseV;
newPTvert(occludedI,:) = [];
faceVInd = 1:1:size(baseV,1);
faceVInd(occludedI) = 0;

% remove all the faces that inpolygon vertices make
% keep track of these faces
TRrem= [];
for i= 1:length(occludedI)
    TRrem = [TRrem;baseTR.ConnectivityList(baseVA{occludedI(i)},:)];
end
TRrem = unique(TRrem,'rows');
baseF = unique(baseF,'rows','stable');
[common, ind] = intersect(baseF,TRrem,'rows');
newTR = baseF;
newTR(ind,:) = [];

newfaceVInd = zeros(1,length(faceVInd));
for i = 1:length(faceVInd)
    newfaceVInd(i) = nnz(faceVInd(1:i));
end
newfaceVInd(1) = 1;
newIndices = reIndexing(newTR,newfaceVInd);
figure;
trimesh(newIndices,newPTvert(:,1),newPTvert(:,2),newPTvert(:,3),...
				'FaceVertexCData',[1,1,1],'EdgeColor',[0,0,0],'FaceColor','flat',...
            'FaceAlpha',0)
%% (4) Stitching
%Delanuay triangle
% the boundary of a hole
tbInd = unique(TRrem(:));
tbInd = tbInd(~ismember(tbInd,occludedI));
newtbInd = reIndexing(tbInd,newfaceVInd);

seamV = [baseV(tbInd,:);partV(partB,:)];
seamI = 1:1:length(seamV);
targetVN = baseVN(centerOccluded,:);
w = null(targetVN);
seamprojV = seamV*w;

% Boundary constraint
C = [];

% base boundary
v = 1:1:length(tbInd);
for i = 1:length(newtbInd)
    b = findRows(newIndices,newtbInd(i));
    for j = 1:size(b,1)
        bb = b(j,:);
        edge = bb(ismember(bb,newtbInd));
        if length(edge)==2
            C = [C;[v(newtbInd==edge(1)),v(newtbInd==edge(2))]];
        else
            continue
        end
    end
end

% partboundary
v = length(tbInd)+1:1:length(seamV);
for i = 1:length(partB)
    b = findRows(partIndices,partB(i));
    for j = 1:size(b,1)
        bb = b(j,:);
        edge = bb(ismember(bb,partB));
        if length(edge) ==2
            C = [C;v(partB == edge(1)),v(partB ==edge(2))];
        else
            continue
        end
    end
end
TRI = delaunayTriangulation(seamprojV(:,1),seamprojV(:,2),C);

figure;
trimesh(TRI.ConnectivityList,seamV(:,1),seamV(:,2),seamV(:,3),...
    'FaceVertexCData',[1,1,1],'EdgeColor',[1,0,1],'FaceColor','flat');
hold on
for i = 1:length(C)
plot3((seamV(C(i,:),1))',(seamV(C(i,:),2))',(seamV(C(i,:),3))','b','LineWidth',1);hold on
end

IO = isInterior(TRI);
newTRI = TRI.ConnectivityList(IO,:);
figure;
subplot(1,3,1);
trimesh(newTRI,seamV(:,1),seamV(:,2),seamV(:,3),...
    'FaceVertexCData',[1,1,1],'EdgeColor',[1,0,1],'FaceColor','flat');
hold on
scatter3(partV(partB,1),partV(partB,2),partV(partB,3),'r','filled');
scatter3(baseV(tbInd,1),baseV(tbInd,2),baseV(tbInd,3),'b','filled');


% order the faces so that they are in the same direction as centerVN
% (facing front of a face) 
ordered_newTRI = zeros(size(newTRI));
for i = 1:length(newTRI)
    if newTRI(i,1)<=length(tbInd) % part of baseMesh
        p1 = baseV(tbInd(newTRI(i,1)),:);
    else % part of partMesh
        p1 = partV(partB(newTRI(i,1)-length(tbInd)),:);
    end
    
    if newTRI(i,2)<=length(tbInd) % part of baseMesh
        p2 = baseV(tbInd(newTRI(i,2)),:);
    else % part of partMesh
        p2 = partV(partB(newTRI(i,2)-length(tbInd)),:);
    end
    
    if newTRI(i,3)<=length(tbInd) % part of baseMesh
        p3= baseV(tbInd(newTRI(i,3)),:);
    else % part of partMesh
        p3 = partV(partB(newTRI(i,3)-length(tbInd)),:);
    end
    v1 = p2-p1;
    v2 = p3-p1;
    
    normal = cross(v1,v2);
    % compare to the normal of a nearby face in the basemesh
    direction = dot(targetVN,normal);
    if direction<0
        ordered_newTRI(i,:) = [newTRI(i,1),newTRI(i,3),newTRI(i,2)]; %[p1,p3,p2]; %clockwise
    else
        ordered_newTRI(i,:) = [newTRI(i,1),newTRI(i,2),newTRI(i,3)]; %[p1,p2,p3]; % counterclockwise
    end
end

subplot(1,3,2);
trimesh(newIndices,newPTvert(:,1),newPTvert(:,2),newPTvert(:,3),...
				'FaceVertexCData',[1,1,1],'EdgeColor',[0,0,0],'FaceColor','flat',...
            'FaceAlpha',0);
        hold on
        trimesh(ordered_newTRI,seamV(:,1),seamV(:,2),seamV(:,3),...
    'FaceVertexCData',[1,1,1],'EdgeColor',[1,0,1],'FaceColor','flat');
trimesh(partIndices,partV(:,1),partV(:,2),partV(:,3));

%% (5) Reindexing the stitched mesh vertices/faces
new_partInd = length(newPTvert)+1:1:length(newPTvert)+size(partV,1);
re_partIndices = reIndexing(partIndices,new_partInd);

remeshed_Indices = [newIndices;re_partIndices];
remeshed_PTvertices = [newPTvert;partV];

allind = 1:1:size(remeshed_PTvertices,1);

range = ismember(remeshed_PTvertices(:,1), baseV(tbInd,1))...
    & ismember(remeshed_PTvertices(:,2), baseV(tbInd,2))...
    & ismember(remeshed_PTvertices(:,3), baseV(tbInd,3));

targetBoundaryInd= allind(range);
partBoundaryInd = new_partInd(partB);

new_seamInd = [targetBoundaryInd,partBoundaryInd];
re_newTRI = reIndexing(ordered_newTRI,new_seamInd);
remeshed_Indices = [remeshed_Indices;re_newTRI];

subplot(1,3,3);
trimesh(remeshed_Indices,remeshed_PTvertices(:,1),remeshed_PTvertices(:,2),...
    remeshed_PTvertices(:,3),'FaceVertexCData',[1,1,1],'EdgeColor',[0,0,0],'FaceColor','flat');
hold on
scatter3(remeshed_PTvertices(targetBoundaryInd,1),remeshed_PTvertices(targetBoundaryInd,2),...
    remeshed_PTvertices(targetBoundaryInd,3),'filled','b');

figure;
trimesh(remeshed_Indices,remeshed_PTvertices(:,1),remeshed_PTvertices(:,2),...
    remeshed_PTvertices(:,3),'FaceVertexCData',[1,1,1],'EdgeColor',[0,0,0],'FaceColor','flat',...
            'FaceAlpha',0);
Lighting;
set(gcf,'position',[1 1 1.15 1.15].*get(gcf,'position'));
 
%% (6) Smoothing
smoothed_PTvertices = MeshSmoothing(remeshed_Indices,remeshed_PTvertices,targetBoundaryInd,partBoundaryInd,smoothingArea);

figure;
subplot(1,3,1);
trimesh(remeshed_Indices,smoothed_PTvertices(:,1),smoothed_PTvertices(:,2),...
    smoothed_PTvertices(:,3),'EdgeColor',[0,0,0],'FaceColor',[1,1,1],'FaceAlpha',0);
hold on;
scatter3(targetCoord(1),targetCoord(2),targetCoord(3),'filled','b');
set(gca,'CameraViewAngle',cameraViewAngle,'CameraTarget',cameraTarget,'CameraPosition',cameraPos);
axis equal
view (0,90)
set(gcf,'color','w');
axis off
set(gcf,'position',[1 1 1.15 1.15].*get(gcf,'position'));

subplot(1,3,2);
trimesh(remeshed_Indices,smoothed_PTvertices(:,1),smoothed_PTvertices(:,2),...
    smoothed_PTvertices(:,3),'EdgeColor','none','FaceColor',[0.65,0.65,0.65]);
Lighting;
set(gcf,'position',[1 1 1.15 1.15].*get(gcf,'position'));

subplot(1,3,3);
trimesh(remeshed_Indices,smoothed_PTvertices(:,1),smoothed_PTvertices(:,2),...
    smoothed_PTvertices(:,3),'EdgeColor','none','FaceColor',[0.65,0.65,0.65]); hold on
Lighting;
set(gcf,'position',[1 1 1.15 1.15].*get(gcf,'position'));
camorbit(90,0,'camera');

newMesh{1,t} = remeshed_Indices;
newMesh{2,t} = smoothed_PTvertices;
newMesh{3,t} = new_partInd;
newMesh{4,t} = part_string;
newMesh{5,t} = targetVN;
newMesh{6,t} = targetCoord;
newMesh{7,t} = targetBoundaryInd;
end
