%% Init
radar = 0; % plot radargram as a subplot with the slope plot
line = 0; % 0 - 19_11, 1 - 18_10-11, 2 - 18_6-7-8

%% Load
load('CustomColormap.mat')
if line == 0
    load("19_11_full.mat")
    bed = readmatrix('19_11_bot.csv');
    bed = bed(:,4)';
elseif line == 1
    A = load("18_10_full.mat");
    B = load("18_11_full.mat");
else
    A = load('19_11_hill_slope.mat');
    B = load("19_11_lake_slope.mat");
    C = load("19_11_bump_slope.mat");

end


%% Graphing
if line == 0
    %Slope
    test1 = -slopegrid;
    test1(test1 == 5) = NaN;
    yd = max(Elevation) - Time/2 * 1.69e8;
    xd = data_x*1000;
    %
    sourceSize = size(test1);
    Data_size = size(Data);
    targetSize = [Data_size(1)-window+1,Data_size(2)];
    [X_samples,Y_samples] = meshgrid(linspace(1,sourceSize(2),targetSize(2)), linspace(1,sourceSize(1),targetSize(1)));
    source_resized_to_target_size = interp2(test1, X_samples, Y_samples);
    full_target = [nan((window-1)/2,size(source_resized_to_target_size,2)); source_resized_to_target_size];
    full_target = vertcat(full_target,nan((window-1)/2,size(source_resized_to_target_size,2)));
    figure(1);
    sgtitle('19 _ 11, Wndw - 100, Flow -->>')
    if radar == 1
        ax1 = subplot(211);
    else
        ax1 = axes;
    end
    surf(xd,yd,full_target,'EdgeColor','none','FaceColor','interp')
    view(180,90)
    shading interp
    clim([-1,1])
    colormap(ax1,flipud(CustomColormap))
    cbar = colorbar(ax1);
    set(ax1, 'YDir', 'reverse');
    cbar.Label.String = 'Reflector Slope, Degrees';
    xlabel(ax1,'Distance, km')
    ylabel(ax1,'Elevation from sea level, m')
    hold on
    plot(ax1,xd,Elevation,'k','LineWidth',3)
    plot(ax1,xd, (Elevation - bed(end:-1:1)),'k',"LineWidth",3)
    % Radargram
    if radar == 1
        ax2 = subplot(212);
        imagesc(ax2,xd, yd(end:-1:1), Data)
        colormap(ax2,'bone')
        colorbar
        set(ax2, 'XDir', 'reverse');
        xlabel(ax2,'Distance, km')
        ylabel(ax2,'Elevation from sea level, m')
        linkaxes([ax1, ax2], 'x','y');
    end
end