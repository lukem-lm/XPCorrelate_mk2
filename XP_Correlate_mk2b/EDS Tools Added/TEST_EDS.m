% PLOT_TIFF_CHANNELS Loads a TIFF file, extracts data from all 4 channels, and plots each channel.

% Define the filename and filepath
tiff_filename = ['EDS_Nb.tiff'];
filepath = '/Users/lukemulholland/Documents/MATLAB/XPCorrelate/Results_11_03_Onwards/'; % Replace with your filepath

% Load the TIFF file
tiff_file = fullfile(filepath, tiff_filename);
C = imread(tiff_file, 'tiff');

% Check the number of channels
num_channels = size(C, 3);
disp(['Number of channels in TIFF file: ', num2str(num_channels)]);

% Display statistics for each channel
for i = 1:num_channels
    channel_data = C(:, :, i);
    disp(['Channel ', num2str(i), ' - Min: ', num2str(min(channel_data(:))), ...
          ', Max: ', num2str(max(channel_data(:))), ...
          ', Mean: ', num2str(mean(channel_data(:)))]);
end

% Plot each channel
figure;
for i = 1:num_channels
    subplot(2, 2, i);
    imagesc(C(:, :, i));
    title(['Channel ', num2str(i)]);
    colorbar;
    axis image;
end

sgtitle(['TIFF File: ', tiff_filename], 'Interpreter', 'none');

