% Load the saved results file
load("Pedro_Reformatted_05.03.25._XPCorrelate_results_Many_Points21-Mar-2025.mat");

% Define the file paths and names for the EDS maps
eds_files = {'EDS_Cr.tiff', 'EDS_Mo.tiff', 'EDS_V.tiff', 'EDS_Nb.tiff'}; % Cr first for reference
filepath = '/Users/lukemulholland/Documents/MATLAB/XPCorrelate/Results_25_3_25/';
resultsdir = filepath;

% Define the channels for each EDS map
eds_channels = struct(...
    'EDS_Mo', 2, ...       % Mo: Channel 2
    'EDS_V', [2, 3], ...   % V: Channels 2 and 3
    'EDS_Nb', [1, 3], ...  % Nb: Channels 1 and 3
    'EDS_Cr', 1 ...        % Cr: Channel 1
);

% Define bins for declination angle
phi_binEdges = 0:5:90;
phi_binCenters = (phi_binEdges(1:end-1) + phi_binEdges(2:end)) / 2;

%% Get Cr Alignment
fprintf('Processing reference file: EDS_Cr.tiff...\n');
[Cr_C, Cr_XC, Cr_YC] = f_loadEPMA('EDS_Cr.tiff', filepath);
Cr_C = Cr_C(:, :, eds_channels.EDS_Cr);

% Get reference points from Cr alignment
[datastack, ref_points] = Fix_EDS_Distortion(datastack, Cr_C, Cr_XC, Cr_YC, 1, 0);
datastack.EPMA_Cr = datastack.EPMAO; % Save Cr map

%% Process other elements using the same reference points
for i = 2:length(eds_files) % Start from 2 since Cr is already processed
    eds_name = eds_files{i};
    fprintf('Processing %s using Cr reference points...\n', eds_name);

    % Load and extract channels
    [C, XC, YC] = f_loadEPMA(eds_name, filepath);
    if strcmp(eds_name, 'EDS_Mo.tiff')
        C = C(:, :, eds_channels.EDS_Mo);
    elseif strcmp(eds_name, 'EDS_V.tiff')
        C = mean(C(:, :, eds_channels.EDS_V), 3);
    elseif strcmp(eds_name, 'EDS_Nb.tiff')
        C = mean(C(:, :, eds_channels.EDS_Nb), 3);
    end

    % Apply the same transformation as Cr
    [locepmacorr(:, 1), locepmacorr(:, 2)] = transformPointsForward(ref_points.tform, datastack.X(:), datastack.Y(:));
    
    % Interpolation
    Cinterp = griddedInterpolant(XC, YC, C, 'linear');
    Onew = Cinterp(locepmacorr(:, 1), locepmacorr(:, 2));
    OnewG = f_gridify_vector(Onew, size(datastack.X, 1), size(datastack.Y, 2))';
    
    % Save aligned map
    switch eds_name
        case 'EDS_Mo.tiff'
            datastack.EPMA_Mo = OnewG;
        case 'EDS_V.tiff'
            datastack.EPMA_V = OnewG;
        case 'EDS_Nb.tiff'
            datastack.EPMA_Nb = OnewG;
    end
end

%% Generate plots for all elements (both hardness and modulus)
for i = 1:length(eds_files)
    generate_combined_plots(datastack, eds_files{i}, phi_binEdges, phi_binCenters, resultsdir);
end

%% Main alignment 
function [datastack, ref_points] = Fix_EDS_Distortion(datastack, C, XC, YC, cropq, EPMArefq)
    % Plot hardness map for reference selection
    figure('Position', [100, 100, 800, 600]);
    contourf(datastack.X, datastack.Y, datastack.H, 455, 'LineColor', 'None');
    axis image;
    colormap(jet);
    colorbar;
    title('Select Reference Points in Hardness Map (4+ points)');
    [x_og, y_og] = getpts;
    close;
    
    % Plot EPMA data for reference selection
    figure('Position', [100, 100, 800, 600]);
    if ndims(C) == 2
        C_plot = C;
    else
        C_plot = C(:, :, 1);
    end
    contourf(XC, YC, C_plot, 45, 'LineColor', 'None');
    axis image;
    colormap(jet);
    colorbar;
    title('Select Corresponding Points in EPMA Map');
    [x_transREF, y_transREF] = getpts;
    close;
    
    % Store reference points
    ref_points.movingPoints = [x_og, y_og];
    ref_points.fixedPoints = [x_transREF, y_transREF];
    ref_points.tform = fitgeotrans(ref_points.movingPoints, ref_points.fixedPoints, 'affine');
    
    % Transformation
    [locepmacorr(:, 1), locepmacorr(:, 2)] = transformPointsForward(ref_points.tform, datastack.X(:), datastack.Y(:));

    % Interpolation
    Cinterp = griddedInterpolant(XC, YC, C, 'linear');
    Onew = Cinterp(locepmacorr(:, 1), locepmacorr(:, 2));
    datastack.EPMAO = f_gridify_vector(Onew, size(datastack.X, 1), size(datastack.Y, 2))'; 
end

%% Plot generation function 
function generate_combined_plots(datastack, eds_name, phi_binEdges, phi_binCenters, resultsdir)
    % Extract element data
    switch eds_name
        case 'EDS_Mo.tiff'
            eds_data = datastack.EPMA_Mo;
            element_name = 'Mo';
        case 'EDS_V.tiff'
            eds_data = datastack.EPMA_V;
            element_name = 'V';
        case 'EDS_Nb.tiff'
            eds_data = datastack.EPMA_Nb;
            element_name = 'Nb';
        case 'EDS_Cr.tiff'
            eds_data = datastack.EPMA_Cr;
            element_name = 'Cr';
    end

    %% 1. Registered EDS map
    fig1 = figure('Name', [element_name ' EDS Map'], 'Position', [100, 100, 800, 600]);
    contourf(datastack.X, datastack.Y, eds_data, 45, 'LineColor', 'None');
    colormap(jet);
    colorbar;
    title([element_name ' Registered EDS Map']);
    xlabel('\mum'); ylabel('\mum');
    axis image;
    save_figure(fig1, element_name, 'Registered_Map', resultsdir);

    %% 2. Scatter plots 
    % Hardness scatter plot
    fig2a = figure('Name', [element_name ' vs Hardness'], 'Position', [100, 100, 800, 600]);
    scatter(eds_data(:), datastack.H(:), 20, 'filled', 'MarkerFaceAlpha', 0.3);
    xlabel([element_name ' Concentration (wt.%)']);
    ylabel('Hardness (GPa)');
    title([element_name ' vs Hardness']);
    grid on;
    save_figure(fig2a, element_name, 'Scatter_Hardness', resultsdir);
    
    % Modulus scatter plot
    fig2b = figure('Name', [element_name ' vs Modulus'], 'Position', [100, 100, 800, 600]);
    scatter(eds_data(:), datastack.M(:), 20, 'filled', 'MarkerFaceAlpha', 0.3);
    xlabel([element_name ' Concentration (wt.%)']);
    ylabel('Modulus (GPa)');
    title([element_name ' vs Modulus']);
    grid on;
    save_figure(fig2b, element_name, 'Scatter_Modulus', resultsdir);

    %% 3. Discrete analysis
    numBins = 10;
    binEdges = linspace(min(eds_data(:)), max(eds_data(:)), numBins + 1);
    binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
    
    % Calculate statistics for both hardness and modulus
    [meanHardness, stdHardness, meanModulus, stdModulus] = deal(zeros(numBins, 1));
    for i = 1:numBins
        binIdx = eds_data(:) >= binEdges(i) & eds_data(:) < binEdges(i+1);
        meanHardness(i) = mean(datastack.H(binIdx));
        stdHardness(i) = std(datastack.H(binIdx));
        meanModulus(i) = mean(datastack.M(binIdx));
        stdModulus(i) = std(datastack.M(binIdx));
    end
    
    % Create results tables
    hardnessTable = table(binEdges(1:end-1)', binEdges(2:end)', meanHardness, stdHardness, ...
        'VariableNames', {'Conc_Lower', 'Conc_Upper', 'Mean_Hardness_GPa', 'Std_Hardness_GPa'});
    modulusTable = table(binEdges(1:end-1)', binEdges(2:end)', meanModulus, stdModulus, ...
        'VariableNames', {'Conc_Lower', 'Conc_Upper', 'Mean_Modulus_GPa', 'Std_Modulus_GPa'});
    
    disp([element_name ' Binned Hardness Results:']);
    disp(hardnessTable);
    disp([element_name ' Binned Modulus Results:']);
    disp(modulusTable);

    % Plot binned data - Hardness
    fig3a = figure('Name', [element_name ' Binned Hardness'], 'Position', [100, 100, 800, 600]);
    bar(binCenters, meanHardness, 'FaceColor', [0.7 0.7 0.7]);
    hold on;
    errorbar(binCenters, meanHardness, stdHardness, 'k.', 'LineWidth', 1.5);
    xlabel([element_name ' Concentration (wt.%)']);
    ylabel('Hardness (GPa)');
    title([element_name ' Binned Hardness Analysis']);
    grid on;
    save_figure(fig3a, element_name, 'Binned_Hardness', resultsdir);
    
    % Plot binned data - Modulus
    fig3b = figure('Name', [element_name ' Binned Modulus'], 'Position', [100, 100, 800, 600]);
    bar(binCenters, meanModulus, 'FaceColor', [0.7 0.7 0.7]);
    hold on;
    errorbar(binCenters, meanModulus, stdModulus, 'k.', 'LineWidth', 1.5);
    xlabel([element_name ' Concentration (wt.%)']);
    ylabel('Modulus (GPa)');
    title([element_name ' Binned Modulus Analysis']);
    grid on;
    save_figure(fig3b, element_name, 'Binned_Modulus', resultsdir);

    %% 4. 3D plots
    % Hardness 3D plot
    meanHardness3D = zeros(length(phi_binCenters), numBins);
    for i = 1:length(phi_binCenters)
        for j = 1:numBins
            idx = eds_data(:) >= binEdges(j) & eds_data(:) < binEdges(j+1) & ...
                  datastack.Phi(:)*180/pi >= phi_binEdges(i) & ...
                  datastack.Phi(:)*180/pi < phi_binEdges(i+1);
            if sum(idx) > 5
                meanHardness3D(i,j) = mean(datastack.H(idx));
            end
        end
    end
    
    fig4a = figure('Name', [element_name ' 3D Hardness'], 'Position', [100, 100, 1000, 800]);
    surf(binCenters, phi_binCenters, smoothdata(meanHardness3D,2,'movmean',3), 'EdgeColor','none');
    colormap(jet);
    shading interp;
    xlabel([element_name ' Concentration (wt.%)']);
    ylabel('Declination Angle (deg)');
    zlabel('Hardness (GPa)');
    title([element_name ' vs Angle vs Hardness']);
    colorbar;
    view(3);
    grid on;
    save_figure(fig4a, element_name, '3D_Hardness', resultsdir);
    
    % Modulus 3D plot
    meanModulus3D = zeros(length(phi_binCenters), numBins);
    for i = 1:length(phi_binCenters)
        for j = 1:numBins
            idx = eds_data(:) >= binEdges(j) & eds_data(:) < binEdges(j+1) & ...
                  datastack.Phi(:)*180/pi >= phi_binEdges(i) & ...
                  datastack.Phi(:)*180/pi < phi_binEdges(i+1);
            if sum(idx) > 5
                meanModulus3D(i,j) = mean(datastack.M(idx));
            end
        end
    end
    
    fig4b = figure('Name', [element_name ' 3D Modulus'], 'Position', [100, 100, 1000, 800]);
    surf(binCenters, phi_binCenters, smoothdata(meanModulus3D,2,'movmean',3), 'EdgeColor','none');
    colormap(jet);
    shading interp;
    xlabel([element_name ' Concentration (wt.%)']);
    ylabel('Declination Angle (deg)');
    zlabel('Modulus (GPa)');
    title([element_name ' vs Angle vs Modulus']);
    colorbar;
    view(3);
    grid on;
    save_figure(fig4b, element_name, '3D_Modulus', resultsdir);
end

%% Save figures
function save_figure(figHandle, element_name, plot_type, resultsdir)
    if ~exist(resultsdir, 'dir')
        mkdir(resultsdir);
    end
    
    filename_png = fullfile(resultsdir, [element_name '_' plot_type '.png']);
    filename_fig = fullfile(resultsdir, [element_name '_' plot_type '.fig']);
    
    % Save as PNG 
    exportgraphics(figHandle, filename_png, 'Resolution', 300);
    disp(['Saved: ' filename_png]);
    
    % Save as MATLAB figure
    savefig(figHandle, filename_fig);
    disp(['Saved: ' filename_fig]);
    
    close(figHandle); 
end