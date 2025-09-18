%% Plot ALL identity and emotions 
% Construct data matrix for 12 id x 7 emotions

num_faces = 12;
num_vertices = 936;
individuals = {'Elias', 'Neptune', 'Seojin', 'Sophie', 'Dan', 'Sreyas', 'Younah', 'Ashley', 'Josh', 'Kedar', 'Ani', 'Tony'};
emotions = {'neutral', 'happiness_4',  'surprise_4', 'sadness_4', 'disgust_4', 'anger_4','fear_4'};
emotions_scaled ={'neutral', 'happiness_4',  'surprise_4', 'sadness_4', 'disgust_4', 'anger_4','fear_4'};

all_faces_data = zeros(num_faces * length(emotions_scaled), num_vertices * 3);

% Map for assigning distinct colors for identities
identity_colors = lines(num_faces);  % Use distinct colors for each identity

counter = 1;
for i = 1:num_faces
    for j = 1:length(emotions_scaled)
        face = scaledCoordsMax_all_updated.([individuals{i} '_' emotions_scaled{j}]);
        all_faces_data(counter, :) = reshape(face, 1, []);  % Flatten and store
        counter = counter + 1;
    end
end
%% Apply PCA

% Mean subtraction (centering the data)
mean_face = mean(all_faces_data, 1);
centered_faces = all_faces_data - mean_face;

% Normalize the data (optional but recommended for PCA)
standardized_faces = centered_faces ./ std(centered_faces);

% Apply PCA
[coeff, score, latent] = pca(standardized_faces);

% Explained variance for the combined set
explained_variance = latent / sum(latent);
disp('Explained variance ratio for each PC:');
disp(explained_variance);

% % Create VideoWriter object to save video as an .mp4 file
% videoFile = VideoWriter('3D_PCA_Rotation_ALL2.mp4', 'MPEG-4');
% videoFile.FrameRate = 10;  % Set frame rate (frames per second)
% open(videoFile);

% 3D Scatter plot with varying color intensity for emotions (with legend)
figure;
hold on;
legend_handles = [];  % Initialize an array for legend handles

for i = 1:num_faces
    base_color = identity_colors(i, :);  % Neutral face color for each identity
    
    % Create a plot for the neutral face only for legend purposes
    h = scatter3(NaN, NaN, NaN, 100, 'filled', 'MarkerFaceColor', base_color);  % Dummy plot for the legend
    legend_handles = [legend_handles h];  % Store the handle for the legend
    
    for j = 1:length(emotions)
        % Cap intensity at 0.8 to avoid excessive blending
        emotion_intensity = 0.9 * (j - 1) / (length(emotions) - 1);  % Scale intensity from 0 to 0.8
        current_color = (1 - emotion_intensity) * base_color + emotion_intensity * [0.9, 0.9, 0.9];  % Blend with gray
        
        % Plot each emotion's projection onto the first three PCs
        index = (i - 1) * length(emotions) + j;
        scatter3(score(index, 1), score(index, 2), score(index, 3), 90, 'filled', 'MarkerFaceColor', current_color);
        
        % Label only the neutral face with the identity name
        if strcmp(emotions{j}, 'neutral')
            text(score(index, 1), score(index, 2), score(index, 3), individuals{i}, 'FontSize', 17);
        end
    end
end

xlabel('PC 1');
ylabel('PC 2');
zlabel('PC 3');

view(3); % Default 3D view
title('3D PCA of Identity and Emotions', FontSize=17);
grid on;
centroid = mean(score(:, 1:3));  % Mean of all data points in the first 3 PCs

% Add the legend with identity names
all_legend_handles = [legend_handles];
all_legend_labels = individuals;
legend(all_legend_handles, all_legend_labels, 'Location', 'bestoutside');  % or 'best'

% % Rotate and capture frames
% for angle = 0:360  % Rotate from 0 to 360 degrees
%     view(angle, 30);  % Adjust the view angle (azimuth = angle, elevation = 30)
%     frame = getframe(gcf);  % Capture the current figure as a frame
%     writeVideo(videoFile, frame);  % Write the frame to the video file
% end
% 
% % Close the video file
% close(videoFile);
% disp('3D plot saved as video.');
%% Identity & Emotion Separation Vectors via Covariance Eigenvectors 

num_faces = 12;
% emotions = {'neutral', 'happiness_4', 'sadness_4', 'disgust_4', 'fear_4', 'anger_4', 'surprise_4'};

% 1. Identity separation vector (variation across identities)
identity_centroids = zeros(num_faces, 3);
for i = 1:num_faces
    indices = (i-1)*length(emotions) + 1 : i*length(emotions);  % All emotions for identity i
    identity_centroids(i, :) = mean(score(indices, 1:3), 1);
end

cov_identity = cov(identity_centroids);
[evecs_id, evals_id] = eig(cov_identity);
[~, idx_id] = sort(diag(evals_id), 'descend');
identity_separation_vector = evecs_id(:, idx_id(1));  % Principal direction
identity_separation_vector = identity_separation_vector / norm(identity_separation_vector);

% 2. Emotion separation vector (variation across emotions)
emotion_centroids = zeros(length(emotions), 3);
for j = 1:length(emotions)
    indices = j:length(emotions):num_faces * length(emotions);  % All identities for emotion j
    emotion_centroids(j, :) = mean(score(indices, 1:3), 1);
end

cov_emotion = cov(emotion_centroids);
[evecs_em, evals_em] = eig(cov_emotion);
[~, idx_em] = sort(diag(evals_em), 'descend');
emotion_separation_vector = evecs_em(:, idx_em(1));
emotion_separation_vector = emotion_separation_vector / norm(emotion_separation_vector);

% === Angle Between Separation Vectors ===
dot_product = dot(identity_separation_vector, emotion_separation_vector);
theta_deg = rad2deg(acos(dot_product));
disp(['Angle between identity and emotion vectors: ', num2str(theta_deg, '%.2f'), ' degrees']);
%% 
num_faces = 12;
num_emotions = length(emotions);  % 7

% Use full PCA scores instead of just top 3
% Assume score_full is (84 x N), where N is number of retained PCs
identity_centroids = zeros(num_faces, size(score, 2));
for i = 1:num_faces
    indices = (i-1)*num_emotions + 1 : i*num_emotions;
    identity_centroids(i, :) = mean(score(indices, :), 1);
end

% Identity separation vector via PCA on centroids
cov_identity = cov(identity_centroids);
[evecs_id, evals_id] = eig(cov_identity);
[~, idx_id] = sort(diag(evals_id), 'descend');
identity_vector_full = evecs_id(:, idx_id(1));
identity_vector_full = identity_vector_full / norm(identity_vector_full);

% Emotion separation vector
emotion_centroids = zeros(num_emotions, size(score, 2));
for j = 1:num_emotions
    indices = j:num_emotions:num_faces * num_emotions;
    emotion_centroids(j, :) = mean(score(indices, :), 1);
end

cov_emotion = cov(emotion_centroids);
[evecs_em, evals_em] = eig(cov_emotion);
[~, idx_em] = sort(diag(evals_em), 'descend');
emotion_vector_full = evecs_em(:, idx_em(1));
emotion_vector_full = emotion_vector_full / norm(emotion_vector_full);

% Compute angle between full vectors
dot_product = dot(identity_vector_full, emotion_vector_full);
angle_rad = acos(dot_product);
angle_deg = rad2deg(angle_rad);

fprintf('Angle between full identity and emotion vectors: %.2f degrees\n', angle_deg);

%% 84x84 Distance Matrix in PCA Space 

% Extract only the first 3 principal components
score7D = score(:, 1:7);
score_full = score;
% score3D_z = zscore(score3D, 0, 2);  % <-- INSERTED HERE
% Initialize matrix
pca_distance_matrix_7d = zeros(size(score7D, 1));
pca_distance_matrix_full = zeros(size(score_full, 1));

% Compute pairwise Euclidean distances
for i = 1:size(score7D, 1)
    for j = 1:size(score7D, 1)
        pca_distance_matrix_7d(i, j) = norm(score7D(i, :) - score7D(j, :));
    end
end
for i = 1:size(score_full, 1)
    for j = 1:size(score_full, 1)
        pca_distance_matrix_full(i, j) = norm(score_full(i, :) - score_full(j, :));
    end
end

% Visualize the matrix
identities = {'Elias', 'Neptune', 'Seojin', 'Sophie', 'Dan', 'Sreyas', ...
              'Younah', 'Ashley', 'Josh', 'Kedar', 'Ani', 'Tony'};
tick_positions = 3.5:7:84;
section_bounds = 1:7:85;  % includes the final edge

figure;
imagesc(pca_distance_matrix_full);
colorbar;
title('84x84 Distance Matrix in Full PCA Space');
xlabel('Face index');
ylabel('Face index');

% Set axis ticks and labels
set(gca, 'XTick', tick_positions, 'XTickLabel', identities, ...
         'YTick', tick_positions, 'YTickLabel', identities);
xtickangle(90);
ytickangle(0);
axis square;

% Overlay lines at start/end of each identity group
hold on;
for k = 1:length(section_bounds)
    x = section_bounds(k) - 0.5;  % align with image pixel edge
    y = section_bounds(k) - 0.5;
    xline(x, 'k-', 'LineWidth', 0.5);  % vertical line
    yline(y, 'k-', 'LineWidth', 0.5);  % horizontal line
end
hold off;
%% 84x84 Distance Matrix in PCA Space 

% Extract only first 3 PCs
score7D = score(:, 1:7);
score_full = score;

% Initialize distance matrices
pca_distance_matrix_7d = zeros(size(score7D, 1));
pca_distance_matrix_full = zeros(size(score_full, 1));

% Compute pairwise Euclidean distances
for i = 1:size(score7D, 1)
    for j = 1:size(score7D, 1)
        pca_distance_matrix_7d(i, j) = norm(score7D(i, :) - score7D(j, :));
    end
end
for i = 1:size(score_full, 1)
    for j = 1:size(score_full, 1)
        pca_distance_matrix_full(i, j) = norm(score_full(i, :) - score_full(j, :));
    end
end

% New grouping by emotion
num_identities = 12;
num_emotions = 7;

new_order = [];
for e = 1:num_emotions
    for i = 1:num_identities
        idx = (i - 1) * num_emotions + e;
        new_order = [new_order; idx];
    end
end

pca_distance_matrix_full_reordered = pca_distance_matrix_full(new_order, new_order);

% Emotion labels
emotions = {'neutral', 'happiness_4', 'surprise_4', 'sadness_4', ...
            'disgust_4', 'anger_4', 'fear_4'};

% Tick positions
tick_positions = 6.5:12:84;
section_bounds = 1:12:85;

% Plot
figure;
imagesc(pca_distance_matrix_full_reordered);
colorbar;
title('84x84 Distance Matrix in Full PCA Space (Grouped by Emotion)');
xlabel('Emotion');
ylabel('Emotion');

set(gca, 'XTick', tick_positions, 'XTickLabel', emotions, ...
         'YTick', tick_positions, 'YTickLabel', emotions);
xtickangle(45);
ytickangle(0);
axis square;

% Overlay lines
hold on;
for k = 1:length(section_bounds)
    x = section_bounds(k) - 0.5;
    y = section_bounds(k) - 0.5;
    xline(x, 'k-', 'LineWidth', 0.5);
    yline(y, 'k-', 'LineWidth', 0.5);
end
hold off;
