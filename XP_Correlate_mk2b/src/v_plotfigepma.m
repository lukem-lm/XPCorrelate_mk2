function v_plotfigepma(datastack, resultsdir, ebsdname, resolution)
    if nargin < 4
        resolution = ['-r' num2str(600)];
    end
    
    % Determine Cr concentration bounds dynamically
    cr_min = min(datastack.EPMAO(:)); 
    cr_max = max(datastack.EPMAO(:)); 
    
    % Plot the registered EDS map with dynamic Cr concentration bounds and jet colormap
    figFC = figure;
    figFC.Name = 'FC';
    hplot = contourf(datastack.X, datastack.Y, datastack.EPMAO, 45, 'LineColor', 'None');
    colormap(jet);
    title('Registered EDS Map (Jet Colormap)')
    xlabel('\mum')
    ylabel('\mum')
    axis image
    c = colorbar;
    c.Label.String = 'Cr Concentration /arb units';
    caxis([cr_min, cr_max]);
    figname = ['FC Figure ' ebsdname(1:(max(size(ebsdname) - 4)))];
    print(fullfile(resultsdir, figname), '-dpng', resolution);
    
    % Plot the EPMA BSE map (if available)
    try
        figepbse = figure;
        figepbse.Name = 'EPMABSE';
        hplot = contourf(datastack.X, datastack.Y, datastack.EPMABSE, 45, 'LineColor', 'None');
        title('Final EPMA BSE Map')
        caxis([95 105]) 
        axis image
        c = colorbar;
        c.Label.String = 'arb units';
        figname = ['EPMABSE Figure ' ebsdname(1:(max(size(ebsdname) - 4)))];
        print(fullfile(resultsdir, figname), '-dpng', resolution);
    end
    
    % Plot Cr concentration vs hardness
    figFCVH = figure;
    figFCVH.Name = 'FC vs H';
    try
        scatter(datastack.EPMAO(:), datastack.Hmat(:));
    catch
        scatter(datastack.EPMAO(:), datastack.H(:));
    end
    xlabel('Cr Concentration /arb units')
    ylabel('Hardness /GPa')
    title({'EDS obtained Cr presence against'; 'measured nanoindentation hardness'})
    xlim([cr_min, cr_max]); 
    figname = ['FC vs H Figure ' ebsdname(1:(max(size(ebsdname) - 4)))];
    print(fullfile(resultsdir, figname), '-dpng', resolution);

    % Plot Cr concentration vs declination angle vs hardness
    figFCVVPH = figure;
    figFCVVPH.Name = 'FC vs phi vs H';
    try
        scatter3(datastack.EPMAO(:), datastack.Phi(:) * 180 / pi, datastack.Hmat(:));
    catch
        scatter3(datastack.EPMAO(:), datastack.Phi(:) * 180 / pi, datastack.H(:));
    end
    xlabel('Cr Concentration /arb units')
    ylabel('Declination angle /^{o}')
    title({'EBSD declination angle against EDS Cr concentration'; 'against measured nanoindentation hardness'})
    ylim([0 90])
    xlim([cr_min, cr_max]); 
    zlabel('Hardness /GPa')
    figname = ['FC vs phi vs H Figure ' ebsdname(1:(max(size(ebsdname) - 4)))];
    print(fullfile(resultsdir, figname), '-dpng', resolution);

    % Plot histogram of hardness vs Cr concentration (discrete intervals)
    figHist = figure;
    figHist.Name = 'Hardness vs Cr Concentration Histogram';

    % Define Cr concentration bins
    numBins = 10; 
    binEdges = linspace(cr_min, cr_max, numBins + 1); 
    binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2; 

    % Initialize arrays to store mean hardness and standard deviation for each bin
    meanHardness = zeros(1, numBins);
    stdHardness = zeros(1, numBins);

    % Calculate mean hardness and standard deviation for each bin
    for i = 1:numBins
        binIndices = datastack.EPMAO(:) >= binEdges(i) & datastack.EPMAO(:) < binEdges(i + 1); % Indices for current bin
        if any(binIndices)
            meanHardness(i) = nanmean(datastack.H(binIndices));
            stdHardness(i) = nanstd(datastack.H(binIndices));
        else
            meanHardness(i) = NaN; 
            stdHardness(i) = NaN;
        end
    end

    % Plot histogram with error bars
    bar(binCenters, meanHardness, 'FaceColor', [0, 0.447, 0.741]);
    hold on;
    errorbar(binCenters, meanHardness, stdHardness, 'k.', 'LineWidth', 1.5); 
    hold off;

    % Adjust axes and labels
    xlabel('Cr Concentration /arb units')
    ylabel('Hardness /GPa')
    xlim([cr_min, cr_max]);
    ylim([0, max(meanHardness + stdHardness, [], 'omitnan') + 1]); 
    title('Hardness vs Cr Concentration (Discrete Intervals)')

    % Add legend with mean and standard deviation for each bin
    legendText = cell(numBins, 1);
    for i = 1:numBins
        legendText{i} = sprintf('%.2f - %.2f: %.2f ± %.2f GPa', ...
                                binEdges(i), binEdges(i + 1), meanHardness(i), stdHardness(i));
    end
    legend(legendText, 'Location', 'eastoutside'); 

    % Save the histogram figure
    figname = ['Hardness vs Cr Concentration Histogram ' ebsdname(1:(max(size(ebsdname) - 4)))];
    print(fullfile(resultsdir, figname), '-dpng', resolution);

    % Create a 3D plot comparing declination angle, Cr concentration, and hardness
    fig3D = figure;
    fig3D.Name = '3D Plot: Hardness vs Cr Concentration vs Declination Angle';

    % Define bins for declination angle
    phi_binEdges = 0:5:90; 
    phi_binCenters = (phi_binEdges(1:end-1) + phi_binEdges(2:end)) / 2; 

    % Matrix to store mean hardness
    meanHardness3D = zeros(length(phi_binCenters), length(binCenters));

    % Calculate mean hardness
    for i = 1:length(phi_binEdges) - 1
        for j = 1:length(binEdges) - 1
            binIndices = datastack.Phi(:) * 180 / pi >= phi_binEdges(i) & ...
                         datastack.Phi(:) * 180 / pi < phi_binEdges(i + 1) & ...
                         datastack.EPMAO(:) >= binEdges(j) & datastack.EPMAO(:) < binEdges(j + 1);
            
            if any(binIndices)
                meanHardness3D(i, j) = nanmean(datastack.H(binIndices));
            else
                meanHardness3D(i, j) = NaN;
            end
        end
    end

    % Create the 3D plot
    surf(binCenters, phi_binCenters, meanHardness3D);
    colormap(jet);
    shading interp; 

    % Add labels and title
    xlabel('Cr Concentration /arb units');
    ylabel('Declination Angle /^{o}');
    zlabel('Hardness /GPa');
    title('3D Plot: Hardness vs Cr Concentration vs Declination Angle');
    colorbar;
    c.Label.String = 'Hardness /GPa';
    view(3); 
    grid on;

    % Save the 3D plot figure
    figname = ['3D Plot Hardness vs Cr vs Phi ' ebsdname(1:(max(size(ebsdname) - 4)))];
    print(fullfile(resultsdir, figname), '-dpng', resolution);

    close all;
end