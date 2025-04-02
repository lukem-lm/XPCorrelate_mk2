%% Clear workspace and add paths
clear all;
close all;
clc;
home;

% Add necessary paths (if any external functions are needed)
addpath(genpath('src'));
addpath(genpath('external'));

%% User inputs
filepath = '/Users/lukemulholland/MATLAB/XPCorrelate/';  % Location of files
epmaname = 'EDS_Cr.tiff'; % Name of EPMA file
epmabsename = 'SE_450X_b.tiff'; % Name of EPMA BSE file (if needed)
cropq = 0; % Crop the EPMA data? 1 = yes, 0 = no
EPMArefq = 0; % Reference mode flag (0 = standard, 1 = BSE correction)

%% Load EPMA data
[C, XC, YC] = f_loadEPMA(epmaname, filepath, epmabsename);

%C = permute(C, [2, 1, 3]); % Transpose C (swap rows and columns for each channel)

% Debug: Check the size of C, XC, and YC
disp('Loaded C dimensions:'); disp(size(C));
disp('Loaded XC dimensions:'); disp(size(XC));
disp('Loaded YC dimensions:'); disp(size(YC));

%% Fix axis orientation 
%XC = XC - min(XC(:)); % Normalize X to start from 0
YC = flipud(YC - min(YC(:))); % Normalize and flip Y-axis

figure;
contourf(C(:, :, 1), 50, 'LineColor', 'none'); % Use Channel 1
title('Raw EPMA Data - Channel 1');
axis image;
colorbar;
colormap('jet');



% Check monotonicity of XC and YC
disp('Checking monotonicity of XC and YC:');
disp(['Is XC sorted? ', num2str(issorted(XC(1, :)))]);
disp(['Is YC sorted? ', num2str(issorted(YC(:, 1)))]);

% Ensure XC and YC are monotonically increasing
if ~issorted(XC(1, :))
    XC = sort(XC, 2); % Sort along rows
end
if ~issorted(YC(:, 1))
    YC = sort(YC, 1); % Sort along columns
end

%% Extract C_plot 
if EPMArefq == 0
    C_plot = C(:, :, 1); % Use correct channel
elseif EPMArefq == 1 && size(C, 3) >= 2
    C_plot = C(:, :, 2); % Use BSE map if available
else
    C_plot = C;
end

%Display matrix sizes
disp('Size of C:'); disp(size(C));
disp('Size of C_plot:'); disp(size(C_plot));

% Ensure XC, YC, and C_plot are compatible
if ~isequal(size(XC), size(YC), size(C_plot))
    error('XC, YC, and C_plot must have the same size.');
end

% Plotting
figure;
contourf(XC, YC, C_plot, 45, 'LineColor', 'none');
colormap("jet");
colorbar;
axis image;
title('Reference Selection in EPMA (>4)');
[x_transREF, y_transREF] = getpts;
close all;

%% Create a dummy datastack for EPMA processing
datastack = struct();
datastack.X = XC; % Use XC as X coordinates
datastack.Y = YC; % Use YC as Y coordinates
datastack.H = C_plot; % Use C_plot as hardness data (dummy)

%% Fix EPMA Distortion
[datastack] = f_fixEPMAdistortion(datastack, C, XC, YC, cropq, EPMArefq);

%% Smooth EPMA Data (if needed)
datastack.EPMAOS = smoothdata(datastack.EPMAO, 2, 'gaussian', 5.5);
datastack.EPMAOS = smoothdata(datastack.EPMAOS, 1, 'gaussian', 5.5);
datastack.EPMAO = datastack.EPMAOS;

%% Plot EPMA Data
figure;
contourf(datastack.X, datastack.Y, datastack.EPMAO, 45, 'LineColor', 'none');
colormap("jet");
colorbar;
axis image;
title('Processed EPMA Data')