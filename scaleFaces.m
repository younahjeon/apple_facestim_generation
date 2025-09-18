function [xyz_scaled] = scaleFaces(xyz_coords, partB)
% scale the input face mesh ((936, 3) vertices) to match the maximum
% average radius of overall faces (0.0662) 

target_radius = 0.0662;
target_radius = 0.081186682318868;
% load("partV_orig_before_scale_Elias.mat")
% xyz_coords = partV_orig_before_scale_Elias;
xyz_coordsB = xyz_coords(partB, :); % boundary coords only
centroid = mean(xyz_coords, 1); % centroid
xy_centered_boundary = xyz_coordsB - centroid; % distances between xy centroid and boundary xy coords
maxDistToCentroidE = max(sqrt(sum(xy_centered_boundary.^2, 2))); % average distance to centroid


xyz_centered = xyz_coords - centroid; % (1) center all xyz coords
scaling_factor = target_radius / maxDistToCentroidE; % (2) find the scaling factor for this face
% disp(mean(sqrt(sum((xyz_centered(partB,:) - mean(xyz_centered, 1)).^2, 2))));
xyz_scaled = xyz_centered * scaling_factor; %  (3) scale the centered coords
% disp(mean(sqrt(sum((xyz_scaled(partB,:) - mean(xyz_scaled, 1)).^2, 2))));

end