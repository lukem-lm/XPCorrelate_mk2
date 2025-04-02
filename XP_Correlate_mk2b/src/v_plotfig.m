function v_plotfig(datastack, resultsdir, ebsdname, saveasfigq)
% Plot fig - CMM script to plot out the various graphs from the data.
try 
    resolution = evalin('base', 'resolution');
catch
    resolution = ['-r' num2str(600)];
end

meanH = nanmean(datastack.H(:));
stdH = nanstd(datastack.H(:));
meanM = nanmean(datastack.M(:));
stdM = nanstd(datastack.M(:));

% Plot nanoindentation hardness map
figure;
try
    meanH = nanmean(datastack.Hmat(:));
    stdH = nanstd(datastack.Hmat(:));
    hplot = contourf(datastack.X, datastack.Y, datastack.Hmat, 455, 'LineColor', 'None');
catch
    hplot = contourf(datastack.X, datastack.Y, datastack.H, 455, 'LineColor', 'None');
end
if meanH > stdH
    caxis([meanH - 2 * stdH, meanH + 2 * stdH])
else 
    caxis([min(hplot(:)), meanH + 1 * stdH])
end

title('Nanoindentation Map')
xlabel('\mum')
ylabel('\mum')
axis image
c = colorbar;
c.Label.String = 'Hardness (GPa)';
figname = ['Hardness Figure ' ebsdname(1:(max(size(ebsdname) - 4)))];
print(fullfile(resultsdir, figname), '-dpng', resolution)
if saveasfigq == 1 
    saveas(gcf, fullfile(resultsdir, figname), 'fig') 
end

% Plot declination angle map
figphi = figure;
figphi.Name = 'Phi';
hplot = contourf(datastack.X, datastack.Y, datastack.Phirefl * 180 / pi(), 45, 'LineColor', 'None');
title('Registered EBSD map')
xlabel('\mum')
ylabel('\mum')
axis image
zlim([0 90])
c = colorbar;
c.Label.String = 'Declination angle /^{o}';
figname = ['Phi Figure paper' ebsdname(1:(max(size(ebsdname) - 4)))];
print(fullfile(resultsdir, figname), '-dpng', resolution)
if saveasfigq == 1 
    saveas(gcf, fullfile(resultsdir, figname), 'fig') 
end

% Plot band contrast map
figbc = figure;
figbc.Name = 'BC';
hplot = contourf(datastack.X, datastack.Y, datastack.BCebsd, 45, 'LineColor', 'None');
title('Registered EBSD map - BC')
xlabel('\mum')
ylabel('\mum')
axis image
c = colorbar;
c.Label.String = 'BC /arb units';
figname = ['BC Figure paper' ebsdname(1:(max(size(ebsdname) - 4)))];
print(fullfile(resultsdir, figname), '-dpng', resolution)
if saveasfigq == 1 
    saveas(gcf, fullfile(resultsdir, figname), 'fig') 
end

% Plot declination angle vs hardness scatter plot
figphiVH = figure;
figphiVH.Name = 'Phi vs H';
scatter(datastack.Phirefl(:) * 180 / pi, datastack.H(:), 'x')
title('Declination angle against measured nanoindentation hardness')
xlabel('Declination angle /^{o}')
ylabel('Hardness /GPa')
ylim([nanmean(datastack.H(:)) - 5 * nanstd(datastack.H(:)), nanmean(datastack.H(:)) + 5 * nanstd(datastack.H(:))])
xlim([0 90])

% Calculate correlation coefficient
phi = datastack.Phirefl(:) * 180 / pi;
H = datastack.H(:);
validIdx = ~isnan(phi) & ~isnan(H);
phi = phi(validIdx);
H = H(validIdx);
corrCoeff = corrcoef(phi, H);
corrText = sprintf('Correlation: %.3f', corrCoeff(1, 2));

% Add correlation coefficient legend
legend(corrText, 'Location', 'best', 'FontSize', 12, 'Color', 'none', 'EdgeColor', 'none');

figname = ['Phi vs H Figure ' ebsdname(1:(max(size(ebsdname) - 4)))];
print(fullfile(resultsdir, figname), '-dpng', resolution)
if saveasfigq == 1 
    saveas(gcf, fullfile(resultsdir, figname), 'fig') 
end

% Plot histogram of hardness vs declination angle
fighist = figure;
fighist.Name = 'Hardness vs Declination Angle Histogram';

% Group declination angles into 5-degree sections
binEdges = 0:5:90; % Bins from 0 to 90 degrees in steps of 5 degrees
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

% Arrays to store mean hardness and standard deviation for each bin
meanHardness = zeros(1, length(binCenters));
stdHardness = zeros(1, length(binCenters));

% Calculate mean hardness and standard deviation
for i = 1:length(binEdges) - 1
    binIndices = phi >= binEdges(i) & phi < binEdges(i + 1); 
    meanHardness(i) = nanmean(H(binIndices)); 
    stdHardness(i) = nanstd(H(binIndices)); 
end

% Plot histogram with error bars
bar(binCenters, meanHardness, 'FaceColor', [0, 0.447, 0.741]); 
hold on;
errorbar(binCenters, meanHardness, stdHardness, 'k.', 'LineWidth', 1.5); 
hold off;

% Adjust axes and labels
xlabel('Declination Angle /^{o}')
ylabel('Hardness /GPa')
xlim([0 90])
ylim([0, max(meanHardness + stdHardness) + 1]); 
title('Hardness vs Declination Angle (5-degree bins)')

% Add legend with mean and standard deviation for each bin
legendText = cell(length(binCenters), 1);
for i = 1:length(binCenters)
    legendText{i} = sprintf('%.0f°: %.2f ± %.2f GPa', binCenters(i), meanHardness(i), stdHardness(i));
end
legend(legendText, 'Location', 'eastoutside'); 

% Save the histogram figure
figname = ['Hardness vs Declination Angle Histogram ' ebsdname(1:(max(size(ebsdname) - 4)))];
print(fullfile(resultsdir, figname), '-dpng', resolution)
if saveasfigq == 1 
    saveas(gcf, fullfile(resultsdir, figname), 'fig') 
end

close all 
end