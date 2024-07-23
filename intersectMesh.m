
function [centerOccluded,occludedI] = intersectMesh()
% intersects the part mesh with the base mesh

global baseF baseV partB baseVA targetCoord baseTR partV baseVN

if ~isempty(targetCoord)
    centerOccluded = nearestNeighbor(baseTR,targetCoord);
else
    centerOccluded = nan;
end

shadowI = nearestNeighbor(baseTR,partV);
shadowF = [];
for i = 1:length(shadowI)
    shadowF = [shadowF;baseTR.ConnectivityList(baseVA{shadowI(i)},:)];
end
shadowF = unique(shadowF,'rows');
shadowV = baseV(unique(shadowF),:);

%centerOccludedVN = baseVN(centerOccluded,:); 
%w = null(targetVN); % project the vertices/faces that are occluded by the part
w = null([0,0,1]);
%w= null(centerOccludedVN); 
% on a 2-dimensional surface
% in other words, flatten the part of potatohead
% that is going to be occluded by the part to caclculate the shadow
shadowprojV = shadowV*w;

hemiInd = baseV(:,1)<=max(shadowV(:,1)) & baseV(:,1)>=min(shadowV(:,1)) &...
    baseV(:,2)<=max(shadowV(:,2)) & baseV(:,2)>=min(shadowV(:,2)) &...
    baseV(:,3)<=max(shadowV(:,3)) & baseV(:,3)>=min(shadowV(:,3));
% Hemisphere of potatohead

hemiV = baseV(hemiInd,:);
hemiprojV = hemiV*w;
inH = inpolygon(hemiprojV(:,1),hemiprojV(:,2),shadowprojV(:,1),shadowprojV(:,2));

in = false(length(baseV),1);
in(hemiInd) = inH;
occludedI = find(in==1);

figure;
subplot(1,2,1);
trimesh(baseF,baseV(:,1),baseV(:,2),baseV(:,3),...
    'FaceVertexCData',[1,1,1],'EdgeColor',[0,0,0],'FaceColor','flat',...
    'FaceAlpha',0); hold on
scatter3(partV(:,1),partV(:,2),partV(:,3),'filled','r');
scatter3(targetCoord(1),targetCoord(2),targetCoord(3),'filled','b','s');
%quiver3(baseV(in,1),baseV(in,2),baseV(in,3),ones(length(occludedI),1)*targetVN(1),ones(length(occludedI),1)*targetVN(2),ones(length(occludedI),1)*targetVN(3),'b');

subplot(1,2,2);
scatter3(partV(partB,1),...
    partV(partB,2),partV(partB,3)); hold on
scatter3(baseV(centerOccluded,1),baseV(centerOccluded,2),baseV(centerOccluded,3),'filled','b');
axis equal
scatter3(baseV(in,1),baseV(in,2),baseV(in,3),'r+');
scatter3(baseV(~in,1),baseV(~in,2),baseV(~in,3),'go');

end


