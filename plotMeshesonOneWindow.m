dotSize = 14;
highlightedDotSize = 20;

individualsToPlot = {'Elias_neutral'};  % add more fields as needed
colors = lines(numel(individualsToPlot));   % auto-colors

figure; hold on
for i = 1:numel(individualsToPlot)
    key = individualsToPlot{i};
    if isfield(scaledCoordsMax_all_updated, key)
        V = scaledCoordsMax_all_updated.(key);  % 936x3
        scatter3(V(:,1), V(:,2), V(:,3), highlightedDotSize, ...
                 'filled', 'MarkerFaceColor', colors(i,:), ...
                 'DisplayName', key);
        % plot3(mean(V(:,1)),mean(V(:,2)),mean(V(:,3)),'x','Color',colors(i,:),'LineWidth',2);
    else
        warning('Missing field: %s', key);
    end
end
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Plot of Selected Individual Meshes');
legend('Location','best','Interpreter','none');
xlim([-0.08 0.08]); ylim([-0.08 0.08]); zlim([-0.08 0.08]);
axis equal vis3d; grid on; rotate3d on; hold off

% Optional export (requires fig2u3d on path)
% fig2u3d(gcf,'save');

%% For the figure
% List of individual mesh field names to plot
individualsToPlot = {'Elias_fear_4', 'Elias_neutral', 'Elias_happiness_4', 'Elias_sadness_4', 'Elias_disgust_4', 'Elias_anger_4', 'Elias_surprise_4', 'Neptune_neutral', 'Seojin_neutral', 'Sophie_neutral', 'Dan_neutral', 'Sreyas_neutral', 'Younah_neutral', 'Ani_neutral', 'Josh_neutral', 'Ashley_neutral', 'Tony_neutral', 'Kedar_neutral'};  % Example
outputDir = './plots_for_paper/';  % <-- Replace with your actual folder path

% Create the directory if it doesn't exist
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

highlightedDotSize = 20;

for i = 1:length(individualsToPlot)
    faceKey = individualsToPlot{i};

    if isfield(scaledCoordsMax_all_updated, faceKey)
        meshCoords = scaledCoordsMax_all_updated.(faceKey);

        % Create figure
        fig = figure('Visible', 'off');
        scatter3(meshCoords(:,1), meshCoords(:,2), meshCoords(:,3), ...
            highlightedDotSize, 'k', 'filled');

        xlim([-0.08, 0.08]);
        ylim([-0.08, 0.08]);
        zlim([-0.08, 0.08]);
        axis equal;
        grid off;
        axis off;

        % Set frontal view (camera pointing toward XY plane from positive Z)
        view([0, 0, 1]);

        % Save figure as JPG
        filename = fullfile(outputDir, [faceKey, '.jpg']);
        exportgraphics(fig, filename, 'BackgroundColor', 'white', 'Resolution', 300);

        close(fig);
    else
        warning('Missing field: %s', faceKey);
    end
end