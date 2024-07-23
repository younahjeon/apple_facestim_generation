function [smoothed_vertices] = MeshSmoothing(faces,vertices,targetBoundaryInd,partBoundaryInd,smoothingArea)
smoothed_vertices = vertices;
remeshedTR = triangulation(faces,vertices);
remeshedVA = vertexAttachments(remeshedTR);

% smooth partBoundary

for i = 1:length(partBoundaryInd)
    neighborF = remeshedTR.ConnectivityList(remeshedVA{partBoundaryInd(i)},:);
    neighborV = unique(neighborF);
    neighborV = neighborV(neighborV ~=partBoundaryInd(i));
    smoothed_vertices(partBoundaryInd(i),:) = mean(vertices(neighborV,:));
end


% the closer the vertices are to the target boundary, the more smoothing
% they will get 
if smoothingArea ~= 0
    [dist,I] = min(pdist2(smoothed_vertices,smoothed_vertices(targetBoundaryInd,:)),[],2);
    smooth_dist = zeros(size(vertices,1),1);
    for i = 1:length(targetBoundaryInd)
        group = find(I== i);
        [sorted_groupdist,sorted_groupInd] = sort(dist(group));
        smoothind = round(smoothingArea/100*length(group));
        smooth_dist(group(sorted_groupInd(1:smoothind))) = sorted_groupdist(1:smoothind);
        % only smooth vertices that are within smoothingArea
        smooth_dist(group(sorted_groupInd(smoothind+1:length(group)))) = nan;
    end
    sigma = max(smooth_dist)/3*3; % standard dev of a normal distribution with mean of 0
    % the larger the number, the greater
    % smoothing will apply to vertices within the
    % smoothingArea
    
    weight = normpdf(smooth_dist,0,sigma);
    weight(isnan(weight)) = 0;
    weight = weight/max(weight); % normalize weight
    
    for i = 1:length(vertices)
        neighborF = remeshedTR.ConnectivityList(remeshedVA{i},:);
        neighborV = unique(neighborF);
        neighborV = neighborV(neighborV~=i);
        smoothed_vertices(i,:) = weight(i)/length(neighborV)*sum(smoothed_vertices(neighborV,:))+ (1-weight(i))*smoothed_vertices(i,:);
    end
end
end
