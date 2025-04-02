function [datastack] = f_fixEPMAdistortion(datastack, C, XC, YC, cropq, EPMArefq)
    % Function to correct EPMA distortion by aligning it to the hardness map.
    % Performs cropping, coordinate normalization, and transformation.

    %% Hardness input
    H = datastack.H;
    X = datastack.X;
    Y = datastack.Y;

    % Plot hardness map for reference selection
    figure;
    contourf(X, Y, H, 455, 'LineColor', 'None');
    axis image;
    caxis([nanmean(H(:)) - 1 * nanstd(H(:)), nanmean(H(:)) + 1 * nanstd(H(:))]);
    title('Reference Selection in Hardness (>4)');
    [x_og, y_og] = getpts; 
    hold on; 
    cropq = 0;
    %% Crop EPMA Data (if requested)
    if cropq == 1
        figure;
        if ndims(C) == 2
            C_plot = C;
        else
            C_plot = C(:, :, 1);
        end
        contourf(XC, YC, C_plot, 45, 'LineColor', 'None');
        title('Cropping EPMA map (select corners)');
        axis image;
        [xc_crop, yc_crop] = getpts;
        close all;

        % Crop data 
        C = C(round(min(yc_crop)):round(max(yc_crop)), round(min(xc_crop)):round(max(xc_crop)), :);
        XC = XC(round(min(yc_crop)):round(max(yc_crop)), round(min(xc_crop)):round(max(xc_crop)));
        YC = YC(round(min(yc_crop)):round(max(yc_crop)), round(min(xc_crop)):round(max(xc_crop)));
    end

    %% Ensure XC and YC are monotonically increasing
    if ~issorted(XC(1, :))
        XC = sort(XC, 2); % Sort along rows
    end
    if ~issorted(YC(:, 1))
        YC = sort(YC, 1); % Sort along columns
    end

    %% Align reference points 
    figure;
    if ndims(C) == 2
        C_plot = C(:, :, 1);
    elseif EPMArefq == 0
        C_plot = C(:, :, 1);
    elseif EPMArefq == 1 && size(C, 3) >= 2
        C_plot = C(:, :, 2);
    else
        C_plot = C;
    end

    % Display sizes
    disp('Size of C:'); disp(size(C));
    disp('Size of C_plot:'); disp(size(C_plot));
    % Check dimensions
    if ~isequal(size(XC), size(YC), size(C_plot))
        error('XC, YC, and C_plot must have the same size.');
    end

    % Plot EPMA data
    contourf(XC, YC, C_plot, 45, 'LineColor', 'None');
    colormap("jet");
    colorbar;
    axis image;
    title('Reference Selection in EPMA (>4)');
    [x_transREF, y_transREF] = getpts;
    close all;

    %% Align plots
    movingPoints = [x_og, y_og];
    fixedPoints = [x_transREF, y_transREF];

    tform = fitgeotrans(movingPoints, fixedPoints, 'affine');
    A = tform.T;

    [locepmacorr(:, 1), locepmacorr(:, 2)] = transformPointsForward(tform, X(:), Y(:));

    %% Interpolation
    disp('Checking sizes before interpolation:');
    disp(['Size of XC: ', num2str(size(XC))]);
    disp(['Size of YC: ', num2str(size(YC))]);
    disp(['Size of C: ', num2str(size(C))]);

    if ~isequal(size(XC), size(YC))
        error('XC and YC must have the same size.');
    end

    if ~isequal(size(XC), size(C(:, :, 1)))
        error('C must have the same spatial dimensions as XC and YC.');
    end

    if ~issorted(XC(1, :)) || ~issorted(YC(:, 1))
        error('XC and YC must be monotonically increasing.');
    end

    % Interpolation for EPMA and BSE data
    if ndims(C) == 2  
        Cinterp = griddedInterpolant(XC, YC, C, 'linear');
        Onew = Cinterp(locepmacorr(:, 1), locepmacorr(:, 2));
        OnewG = f_gridify_vector(Onew, size(X, 1), size(Y, 2))'; 
        datastack.EPMAO = OnewG; 

    elseif size(C, 3) >= 2  
        Cinterp = griddedInterpolant(XC, YC, C(:, :, 1), 'linear');
        Onew = Cinterp(locepmacorr(:, 1), locepmacorr(:, 2));
        OnewG = f_gridify_vector(Onew, size(X, 1), size(Y, 2))'; 

        BSEinterp = griddedInterpolant(XC, YC, C(:, :, 2), 'linear');
        BSEnew = BSEinterp(locepmacorr(:, 1), locepmacorr(:, 2));
        BSEnewG = f_gridify_vector(BSEnew, size(X, 1), size(Y, 2))'; 

        datastack.EPMAO = OnewG;  
        datastack.EPMABSE = BSEnewG;  
    end
end