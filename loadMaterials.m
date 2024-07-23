function [baseF, baseV, partF, partV_orig, partB,refInd,baseTR,baseVN,baseVA,baseFN,...
    partIndices, partTR, partVN, partVA, partFN, partRotAxes] = loadMaterials(id_string,part_string,baseMesh,scale,emotion) 


fname = strcat('FaceParts_',id_string,'.mat'); 
load(fname); 

partF = FaceParts.(part_string){1}; 
partV_orig = FaceParts.(part_string){2}; % scale determines the size of the part 

% emotion 
if ~isempty(emotion)
    partV_orig = changeMeshbyEmotion(emotion,partV_orig,partF);
end

partV_orig = partV_orig* scale;

partB = FaceParts.(part_string){3}; % boundaryInd
if strcmp(baseMesh,'PotatoHead')
    load('newPotatoHead.mat');
    baseV = newhead.vertices;
    baseF = newhead.faces;
else
      error('Wrong Input')
end

baseTR = triangulation(baseF,baseV);
baseVN = vertexNormal(baseTR); 
baseVA = vertexAttachments(baseTR); 
baseFN = faceNormal(baseTR); 
    
partIndices = reIndexing(partF); 
partTR = triangulation(partIndices,partV_orig); 
partVN = vertexNormal(partTR);
partVA = vertexAttachments(partTR);
partFN = faceNormal(partTR); 
figure; trimesh(baseF,baseV(:,1),baseV(:,2),baseV(:,3),'FaceVertexCData',[1,1,1],'EdgeColor',[0,0,0],'FaceColor','flat');
hold on
quiver3(baseV(:,1),baseV(:,2),baseV(:,3),baseVN(:,1),baseVN(:,2),baseVN(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
figure; trimesh(partIndices,partV_orig(:,1),partV_orig(:,2),partV_orig(:,3),'FaceVertexCData',[1,1,1],'EdgeColor',[0,0,0],'FaceColor','flat');
hold on

refInd = FaceParts.(part_string){4}; % reference index of the part 
% for eyes, refInd is the top center of an eye
% for nose, refInd is the nosetip and top of the nose 
% for mouth, refInd is the top center of the mouth

% Find the rotational axes of the part by SVD
[U,S,V] = svd(partV_orig);
% Each column of V is a principal component
partRotAxes = V;

scatter3(partV_orig(partB,1),partV_orig(partB,2),partV_orig(partB,3),'filled','b');
%scatter3(partV(refInd(2),1),partV(refInd(2),2),partV(refInd(2),3),'filled','r'); 
scatter3(partV_orig(refInd,1),partV_orig(refInd,2),partV_orig(refInd,3),'filled','r'); 
h1 = plot3([mean(partV_orig(:,1)),V(1,1)],[mean(partV_orig(:,2)),V(2,1)],[mean(partV_orig(:,3)),V(3,1)],'k')
h2 = plot3([mean(partV_orig(:,1)),V(1,2)],[mean(partV_orig(:,2)),V(2,2)],[mean(partV_orig(:,3)),V(3,2)],'g')
h3 = plot3([mean(partV_orig(:,1)),V(1,3)],[mean(partV_orig(:,2)),V(2,3)],[mean(partV_orig(:,3)),V(3,3)],'m')
quiver3(partV_orig(:,1),partV_orig(:,2),partV_orig(:,3),partVN(:,1),partVN(:,2),partVN(:,3)); 
xlabel('x')
ylabel('y')
zlabel('z')
xlim([min(partV_orig(:,1))-0.05, max(partV_orig(:,1))+0.05])
ylim([min(partV_orig(:,2))-0.05, max(partV_orig(:,2))+0.05])
zlim([min(partV_orig(:,3))-0.05, max(partV_orig(:,3))+0.05])
legend([h1,h2,h3],{'PC1','PC2','PC3'})
end