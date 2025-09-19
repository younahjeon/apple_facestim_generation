function [face_new] = morphMeshes (face1, face2, degree, partB)
% Load all face meshes (scaled)
load("scaledCoordsStructAvg.mat")

% Pick two faces
% face1 = 'Elias';
% face2 = 'Neptune';

% Load vertices
face1_v = scaledCoordsStructAvg.(face1);
face2_v = scaledCoordsStructAvg.(face2);

% Calculate diff
delta = face1_v - face2_v;

% Set the degree of morphing (range from 0 to 1, where 0 results in no change and 1 results in complete morph to face 2)
deg = degree;
face_new = face1_v - delta * deg;

% Scale
face_new = scaleFaces(face_new, partB);

end



