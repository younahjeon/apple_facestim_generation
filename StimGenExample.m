%% 
global parsed_data
FaceParts = struct; 
datafile = 'C:\Users\HANJOO LEE\MATLAB\SeojinE\faceData.txt'; % path to the textfile from FaceCaptureX app 
frame_to_use = 1; % which frame to use? you might want to use a frame where the face has a neutral expression
parsed_data = faceData_readLog(datafile);

load('FacePartsInd.mat');
FaceParts_example = load('FaceParts_example.mat');
load('TriangleIndices.mat');
fields= fieldnames(FacePartsInd);
for i = 1:length(fields)
    cell_to_save = cell(1,4);

    % first cell is the triangular faces corresponding to the specific face
    % part
    cell_to_save{1} = findRows(TriangleIndices', FacePartsInd.(fields{i}),'or+');
   
    % second cell is the actual vertex coordinates
    cell_to_save{2}= parsed_data.vertices{frame_to_use}(FacePartsInd.(fields{i}),:);

    % third cell is vertices that compose the boundary of the part
    cell_to_save{3} = FaceParts_example.FaceParts.(fields{i}){3};

    % fourth cell is reference point
    cell_to_save{4} = FaceParts_example.FaceParts.(fields{i}){4};

    FaceParts.(fields{i}) = cell_to_save;
end

id_string = 'SeojinE'; % 
save(strcat('FaceParts_',id_string), 'FaceParts');


%%

% LOAD PART MESH AND BLANK HEAD
% Load FaceParts_example, which is a structure that contains information
% about different face parts of sophie 
% (1) triangulars face indices, (2) coordinates,
% (3) indices of the outline of the part, 
% (4) and reference index/indices
global targetCoord baseF baseV partV partF partB refInd baseTR baseVN baseVA baseFN partIndices partV_orig_before_scale
baseMesh = 'PotatoHead';
id_string = 'Elias'; % you can change this to any id_string you want to read
part_string = 'ALLs'; % all parts but smaller
scale = 0.9;
emotion = 'neutral';
[baseF, baseV, partF, partV_orig, partB, refInd, baseTR,baseVN,baseVA,baseFN,...
            partIndices, partTR, partVN, partVA, partFN, partRotAxes, partV_orig_before_scale] = loadMaterials(id_string,part_string,baseMesh,scale,emotion);

% partV_orig = scaledCoordsMax.("Elias_Neptune");
% partV_orig = changeMeshbyEmotion(emotion,partV_orig,partF);
% xyz_scaled = scaleFaces(partV_orig, partB);
% partV_orig_before_scale = xyz_scaled;
% scaledCoordsMax_Ashley.(emotion) = partV_orig_before_scale;

% displacement = sqrt(sum((scaledCoordsMax.('SeojinE') - partV_orig_before_scale).^2, 2));
% fprintf('Total Displacement: %.4f\n', sum(displacement) );
% loadMaterials produces two figures, one of the blank head and the second
% of the part. The reference indices of the part are colored red in the
% second figure.
%% (Optional) Apply linear regression
parsed_data = faceData_readLog(datafile);
X = parsed_data.blendshape;   % (nf x nb)
X = X';                     % Transpose to frames × blendshapes
X = [X, ones(size(X,1),1)]; % Add dummy column
% Add dummy column
% X(:, end+1) = 1;

Y = parsed_data.vertices; % e.g., reshape vertex coordinates
nf = length(parsed_data.vertices);
nv = size(parsed_data.vertices{1}, 1);

% Reshape each frame's vertices into a row vector
Y = zeros(nf, nv * 3); % If using x,y,z together

for i = 1:nf
    % Each frame: (nv × 3) -> (1 × nv*3)
    Y(i, :) = reshape(parsed_data.vertices{i}', 1, []);
end
ridgeweight = 0.01;
beta = LinRegressBeta(Y, X, ridgeweight);

%% (Optional) Generating morphed face 
% 1. rotate the part so that the reference point has the same vertex normal
% as the target vertex normal 
% emotions = {'neutral', 'happiness_4', 'sadness_4', 'disgust_4', 'fear_4', 'anger_4', 'surprise_4'};
% numEmotions = length(emotions);
% displacements = zeros(21, 1);  % Pre-allocate for the 21 pairwise displacements
% index = 1;  % To keep track of the displacements array
% 
% % Loop through all emotion pairs
% for i = 1:numEmotions
%     for j = i+1:numEmotions
%         % Calculate displacement between two emotions
%         displacement = sqrt(sum((scaledCoordsMax_Elias.(emotions{i}) - scaledCoordsMax_Elias.(emotions{j})).^2, 2));
%         displacements(index) = sum(displacement);  % Store the mean displacement for each pair
%         disp(emotions{i})
%         disp(emotions{j})
%         disp(sum(displacement))
%         index = index + 1;
%     end
% end
% 
% % Calculate the mean and variance of the displacements
% meanDisplacement = mean(displacements);
% varianceDisplacement = var(displacements);
% 
% % Display the results
% meanDisplacement
% varianceDisplacement
% min(displacements)
% max(displacements)
% partV_orig = morphMeshes('Neptune', 'Elias', 0.5, partB); % morph faces
% displacement = sqrt(sum((scaledCoordsMax.("Neptune") - partV_orig).^2, 2)); % distance between neutral and emotion-applied face
% partV_orig = partV_orig * scale;
% total_displacement = sum(displacement);
% fprintf('Total Displacement: %.4f\n', total_displacement);
% Elias_Neptune = partV_orig;
% scaledCoordsMax.('Elias_Neptune') = scaleFaces(partV_orig, partB);
%% 

refVN = partVN(refInd(2),:);
targetVN = [refVN(1) refVN(2) + 0.00001 refVN(3)];
% targetVN = [-0.06806780653	-0.9027324182	0.4245708401]; 

% targetVN = [0, -0.105, 1.00]; 
%targetVN = [0, -0.1175, 1.35]; % neptune 
partV = rotateMesh(partV_orig,'vertexnormal',refVN,targetVN); 

%%
% 2. Define the location on the blank head 
% where the reference point will be attached to. You can additionally
% define the reference point of the part. For example, ALLs contains two
% reference points

global refInd_optim
refInd_optim = refInd(1);

% [targetCoord targetVN]= targetSelection(baseF,baseV);

% or you can provide targetCoord as below
targetCoord = [-0.002   0.058 0.049]; % happiness (frame=7233)
%%
% 2. do mesh surgery 
smoothingArea = 100; 
t = 1; % step of FaceMorph
newMesh = {};
[newMesh,optim_param]= MeshSurgery(t,part_string,smoothingArea,newMesh,id_string);

TR = triangulation(newMesh{1,t},newMesh{2,t}/max(newMesh{2,t},[],'all'));
filename = strcat(id_string, '_', emotion, '.stl');
stlwrite(TR,filename)


%%
% animate the resulting face with arbitrary blendshapes

t_to_animate = 1; % 3 seconds
fps = 10; 
n_frames = t_to_animate * fps; 
writeMovie = 1;

v = VideoWriter('movie_Seojin.mp4','MPEG-4');
open(v)
load('BETA.mat')

nVertices = 1220;
%X = [0.104969700000000;0.104980100000000;0.114777400000000;0.00303261400000000;0.00296590700000000;0.00666228400000000;0.237118100000000;0.238192200000000;-2.18774000000000e-12;-2.18845600000000e-12;0.0560495000000000;0.0560571100000000;0.102377000000000;0.00716774000000000;0;0;0;0;0.685995100000000;0.685998300000000;0.0816218400000000;0.0811506100000000;0.143178100000000;1.53680700000000e-17;0.151617800000000;0.0437056300000000;0.0116943600000000;0.110847700000000;0.110109400000000;0;0;0.0930285600000000;0.0510150900000000;0.792670200000000;0.786838200000000;0.0177960500000000;0.0201759400000000;0.0215139600000000;1.64507500000000e-05;0.0232451400000000;0.0340266800000000;0.00674136400000000;0.287926900000000;0.917157500000000;0.917360100000000;0.212382300000000;0.198032400000000;0.244017900000000;0.245387400000000;0.279677200000000;0.270702900000000;0.000261358700000000]
X = zeros(52,1);
X(9) = 1;
X(10) = 1;
X(length(X)+1) = 0;
delta = X'* BETA; 
delta = reshape(delta,nVertices,3); 
load('FacePartsInd.mat'); 
ALLs_ind = FacePartsInd.ALLs;
delta_sub = delta(ALLs_ind,:);

neutral = newMesh{2};

vertices_orig = neutral;

vertices_target = vertices_orig;
vertices_target(newMesh{3},:) = vertices_orig(newMesh{3},:) + delta_sub;

vertices_inc = (vertices_target - vertices_orig)/ n_frames; 
figure;
%set(gcf, 'renderer', 'zbuffer');
for i = 1:n_frames
    vertices = vertices_orig + (i-1) * vertices_inc;
   
    trimesh(newMesh{1},vertices(:,1),vertices(:,2),vertices(:,3),...
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
    if i ==1
    set(gcf,'position',[1 1 1.15 1.15].*get(gcf,'position'));
    view (0,90)
    end
    title(sprintf('animation time = %2.2f seconds \nframe = %d',t_to_animate,i));
    drawnow 
    pause(0.01)
    if writeMovie
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end
close(v); % crucial! 


%% Compute the magnitude of vertex displacement 
% Assuming partV, partV2 as face parts before and after change
% 
% displacement = sqrt(sum((partV_orig_before_scale_ - partV_orig_before_scale_Sreyas).^2, 2)); % Euclidean distance between corresponding vertices, calculated per row (vertex)
% normalizedDisplacements = (displacement - min(displacement)) / (max(displacement) - min(displacement));
% normalizedDisplacements = (displacement - 0.00026281) / (0.02943009 - 0.00026281);
% 
% znormalizedDisplacements = abs((displacement - mean(displacement)) / std(displacement)); % z-score
% 
% [sortedDisplacements, originalIndices] = sort(displacement, 'ascend');
% 
% total_displacement = sum(normalizedDisplacements);
% mean_displacement = mean(normalizedDisplacements);
% max_displacement = max(normalizedDisplacements);
% 
% fprintf('Total Displacement: %.4f\n', total_displacement);
% fprintf('Mean Displacement: %.4f\n', mean_displacement);
% fprintf('Max Displacement: %.4f\n', max_displacement);

% % Sort the displacements in descending order and get the indices
% [sorted_displacement, sorted_indices] = sort(displacement, 'descend');
% 
% % Get the top 10 displacements and their corresponding indices
% top_10_displacements = sorted_displacement(1:10);
% top_10_indices = sorted_indices(1:10);
% 
% % Print the results
% fprintf('Top 10 Displacements and their Indices:\n');
% for i = 1:10
%     fprintf('Index: %d, Displacement: %.4f\n', top_10_indices(i), top_10_displacements(i));
% end
% %% % Assuming partV_orig_before_scale_X are the matrices for each identity
% identities = {'Elias', 'Dan', 'Sophie', 'Neptune', 'SeojinE', 'Sreyas', 'Younah'};
% 
% % Initialize arrays to store all displacements
% allDisplacements = zeros(936, 1);
% 
% % Compute pairwise displacements and concatenate them
% for i = 1:length(identities)
%     for j = i+1:length(identities)
%         % Get the vertex matrices for the two identities
%         identity1 = eval(['partV_orig_before_scale_' identities{i}]);
%         identity2 = eval(['partV_orig_before_scale_' identities{j}]);
% 
%         % Compute the Euclidean distance (displacement) between corresponding vertices
%         squaredDifferences = (identity1 - identity2).^2;
%         displacements = sqrt(sum(squaredDifferences, 2));
% 
%         % Concatenate the displacements
%         allDisplacements = allDisplacements + displacements;
%     end
% end
% 
% % Find the global minimum and maximum displacements
% minDisplacement = min(allDisplacements);
% maxDisplacement = max(allDisplacements);
% 
% % Display the results
% fprintf('The minimum displacement across all pairs is: %.8f\n', minDisplacement);% 0.00026281
% fprintf('The maximum displacement across all pairs is: %.8f\n', maxDisplacement); % 0.02943009
% 
% %% 
% 
% % Mapping names to their corresponding data variables
% dataMap = containers.Map({'Elias', 'Dan', 'Sophie', 'Neptune', 'SeojinE', 'Sreyas', 'Younah', 'Elias_2024'}, ...
%     {'partV_orig_before_scale_Elias', 'partV_orig_before_scale_Dan', 'partV_orig_before_scale_Sophie', ...
%     'partV_orig_before_scale_Neptune', 'partV_orig_before_scale_SeojinE', 'partV_orig_before_scale_Sreyas', ...
%     'partV_orig_before_scale_Younah', 'partV_orig_before_scale_Elias_2024'});
% 
% % Specify two IDs to compare
% name1 = 'Elias';
% name2 = 'Dan';
% 
% % Retrieve the corresponding data based on the specified names
% data1 = eval(dataMap(name1));
% data2 = eval(dataMap(name2));
% 
% % Compute the Euclidean distance (displacement) between corresponding vertices
% displacement = sqrt(sum((data1 - data2).^2, 2));
% normalizedDisplacements = (displacement - min(displacement)) / (max(displacement) - min(displacement));
% 
% % dotSize = 5; % Adjust this value to make the dots smaller or larger
% % figure;
% % scatter3(data1(:,1), data1(:,2), data1(:,3), dotSize, normalizedDisplacements, 'filled');
% % xlabel('X');
% % ylabel('Y');
% % zlabel('Z');
% % 
% % title(['Heatmap of Displacements: ' name1 ' vs ' name2]);
% % 
% % colorbar;
% % colormap('jet'); 
% % 
% % xlim([min(data1(:,1)) - 0.05, max(data1(:,1)) + 0.05]);
% % ylim([min(data1(:,2)) - 0.05, max(data1(:,2)) + 0.05]);
% % zlim([min(data1(:,3)) - 0.05, max(data1(:,3)) + 0.05]);
% % 
% % hold off;
% 
% % Sort the displacements in ascending order and get the indices
% [sorted_displacement, sorted_indices] = sort(allDisplacements, 'ascend');
% 
% % Get the top 10 displacements and their corresponding indices
% bottom_50_displacements = sorted_displacement(1:50);
% bottom_50_indices = sorted_indices(1:50);
% 
% top_50_displacements = sorted_displacement(end:-1:end-49);
% top_50_indices = sorted_indices(end:-1:end-49);
% 
% % Define the size of the dots (smaller value for smaller dots)
% dotSize = 5; % Adjust this value to make the dots smaller or larger
% 
% % Initialize the figure window
% figure;
% 
% % Plot the face mesh with color based on displacement values
% scatter3(data1(:,1), data1(:,2), data1(:,3), dotSize, allDisplacements, 'filled');
% hold on;
% 
% % Highlight the top 10 displacements with a different color and larger size
% % highlightedDotSize = 30; % Larger size for top 10 points
% % scatter3(data1(top_50_indices, 1), data1(top_50_indices, 2), data1(top_50_indices, 3), highlightedDotSize, 'r', 'filled');
% % scatter3(data1(bottom_50_indices, 1), data1(bottom_50_indices, 2), data1(bottom_50_indices, 3), highlightedDotSize, 'b', 'filled');
% 
% % Set labels
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% 
% % Set the title
% title(['Heatmap of Displacements: ' name1 ' vs ' name2]);
% 
% % Add a colorbar to show the displacement scale
% colorbar;
% colormap('jet'); % You can change 'jet' to another colormap if preferred
% 
% % Adjust axis limits for better visualization
% xlim([min(data1(:,1)) - 0.05, max(data1(:,1)) + 0.05]);
% ylim([min(data1(:,2)) - 0.05, max(data1(:,2)) + 0.05]);
% zlim([min(data1(:,3)) - 0.05, max(data1(:,3)) + 0.05]);
% 
% % Hold off to stop adding to this plot
% hold off;
% %% 
% load('FacePartsInd.mat');
% load('TriangleIndices.mat'); % Face mesh topology
% % load('SomeParsedFace.mat'); % The full 1220-vertex face
% 
% vertices = parsed_data.vertices{frame_to_use}; % 1220 × 3
% included = false(size(vertices,1),1);
% included(FacePartsInd.ALLs) = true;
% 
% scatter3(vertices(:,1), vertices(:,2), vertices(:,3), 10, included, 'filled');
% axis equal; colorbar; title('Vertices in ALLs (936 out of 1220)');
% 
