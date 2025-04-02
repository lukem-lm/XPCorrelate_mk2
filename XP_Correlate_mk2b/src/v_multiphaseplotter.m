function v_multiphaseplotter(datastack, resultsdir, ebsdname)
% Plots various things that are helpful if dealing with multiple phases and
% looking at phase-property relationships. CHECK THE VALUES FOR PHASE
% DISCRIMINATION.

    % Plot the final phase map
    figphas = figure;
    figphas.Name = 'Phase';
    hplot = contourf(datastack.X, datastack.Y, datastack.phase, 45, 'LineColor', 'None');
    colormap(jet); % Varied to jet colour from default
    colorbar;
    title('Final Phase Map')
    axis image
    c = colorbar;
    c.Label.String = 'Phase';
    figname = ['Phase Figure ' ebsdname(1:(max(size(ebsdname) - 4)))];
    saveas(gcf, fullfile(resultsdir, figname), 'png');
    
    % Plot phase vs hardness (scatter plot)
    figphasVH = figure;
    figphasVH.Name = 'phase vs H';
    scatter(datastack.phase(:), datastack.H(:));
    colormap(jet); % Varied to jet colour from default
    colorbar;
    xlabel('Phase')
    ylabel('Hardness /GPa')
    figname = ['phase vs H Figure ' ebsdname(1:(max(size(ebsdname) - 4)))];
    saveas(gcf, fullfile(resultsdir, figname), 'png');

    % Find values of locations where the following conditions are satisfied
    where1 = datastack.phase > 0.9 & datastack.phase < 1.1; % Phase 1
    where2 = datastack.phase > 1.9 & datastack.phase < 2.1; % Phase 2

    % Get means and standard deviations where H is in that phase
    meanH1 = nanmean(datastack.H(where1)); % Mean hardness for Phase 1
    meanH2 = nanmean(datastack.H(where2)); % Mean hardness for Phase 2
    
    stdevH1 = nanstd(datastack.H(where1)); % Standard deviation for Phase 1
    stdevH2 = nanstd(datastack.H(where2)); % Standard deviation for Phase 2
    
    % Create phase vs mean hardness data with error bars
    phasevH = [1, meanH1, stdevH1; 
               2, meanH2, stdevH2];
    
    % Plot phase vs mean hardness as a bar chart with error bars
    figphasVHM = figure;
    figphasVHM.Name = 'Phase vs Mean of H';
    
    % Define bar positions and width
    barPositions = phasevH(:, 1); % Bar positions (x-axis)
    barWidth = 0.5; % Width of the bars
    bar(barPositions, phasevH(:, 2), barWidth, 'FaceColor', [0, 0.447, 0.741]); % Blue bars
    
    hold on;
    errorbar(barPositions, phasevH(:, 2), phasevH(:, 3), 'k.', 'LineWidth', 1.5); % Error bars
    hold off;
    
    % Add legend with hardness values and errors
    legendText = cell(size(phasevH, 1), 1);
    for i = 1:size(phasevH, 1)
        legendText{i} = sprintf('Phase %d: %.2f ± %.2f GPa', ...
                                phasevH(i, 1), phasevH(i, 2), phasevH(i, 3));
    end
    legend(legendText, 'Location', 'northeastoutside'); % Legend outside the plot
    
    % Adjust axes
    xlabel('Phase')
    ylabel('Hardness /GPa')
    
    % Set x-axis limits with padding
    xPadding = 0.5; % Padding on either side of the bars
    xlim([min(barPositions) - xPadding, max(barPositions) + xPadding]); % X-axis limits with padding
    
    % Set y-axis limits
    ylim([0, max(phasevH(:, 2)) + max(phasevH(:, 3)) + 1]); % Y-axis from 0 to max hardness + error + buffer
    
    % Custom X-axis labels
    xlabelnames = {'Phase 1', 'Phase 2'};
    set(gca, 'xtick', barPositions, 'xticklabel', xlabelnames); % Centered labels
    
    % Save the figure
    figname = ['phase vs H Mean Figure ' ebsdname(1:(max(size(ebsdname) - 4)))];
    saveas(gcf, fullfile(resultsdir, figname), 'png');
    
    % Plot hardness vs modulus (scatter plot)
    figphasMVH = figure;
    figphasMVH.Name = 'Hardness vs Modulus';
    scatter(datastack.M(where1), datastack.H(where1), 15, 'filled', 'r'); % Phase 1
    hold on;
    scatter(datastack.M(where2), datastack.H(where2), 15, 'filled', 'b'); % Phase 2
    hold off;
    xlabel('Modulus /GPa')
    ylabel('Hardness /GPa')
    legend('Phase 1', 'Phase 2', 'Location', 'Best');
    figname = ['Hardness vs Modulus Figure ' ebsdname(1:(max(size(ebsdname) - 4)))];
    saveas(gcf, fullfile(resultsdir, figname), 'png');
end