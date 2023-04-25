%% Init
clear all
close all

radar = 0; % plot radargram as a subplot with the slope plot
line = 1; % 0 - 19_11, 1 - 18_10-11, 2 - 18_6-7-8

%% Load
load('CustomColormap.mat')
if line == 0
    load("19_11_full.mat")
    bed = readmatrix('19_11_bot.csv');
    bed = bed(:,4)';
elseif line == 1
    A = load("18_10_full.mat");
    bed_A = readmatrix("18_10_bot.csv");
    bed_A = bed_A(:,4)';
    B = load("18_11_full.mat");
    bed_B = readmatrix("18_11_bot.csv");
    bed_B = bed_B(:,4)';
else
    A = load('19_11_hill_slope.mat');
    B = load("19_11_lake_slope.mat");
    C = load("19_11_bump_slope.mat");

end


%% Graphing
if line == 0
    %Slope
    slp = -slopegrid;
    slp(slp == 5) = NaN;
    yd = max(Elevation) - Time/2 * 1.69e8;
    xd = data_x*1000;
    %
    sourceSize = size(slp);
    Data_size = size(Data);
    targetSize = [Data_size(1)-window+1,Data_size(2)];
    [X_samples,Y_samples] = meshgrid(linspace(1,sourceSize(2),targetSize(2)), linspace(1,sourceSize(1),targetSize(1)));
    source_resized_to_target_size = interp2(slp, X_samples, Y_samples);
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

if line == 1 || line == 2
    %Slope
    if line == 1
        xd = zeros(1);
        yd = [];
        for i = 1:2
            if i == 1
                S = A;
            else 
                S = B;
            end
            %slope
            slp = -S.slopegrid;
            slp(slp == 6) = NaN;
            yd = [yd; max(S.Elevation) - (S.Time/2 * 1.69e8)];
            xd = [xd (S.data_x*1000 + max(xd))];
            %resize
            sourceSize = size(slp);
            Data_size = size(S.Data);
            targetSize = [Data_size(1)-S.window+1,Data_size(2)];
            [X_samples,Y_samples] = meshgrid(linspace(1,sourceSize(2),targetSize(2)), linspace(1,sourceSize(1),targetSize(1)));
            source_resized_to_target_size = interp2(slp, X_samples, Y_samples);
            full_target = [nan((S.window-1)/2,size(source_resized_to_target_size,2)); source_resized_to_target_size];
            full_target = vertcat(full_target,nan((S.window-1)/2,size(source_resized_to_target_size,2)));
            if i ==1
                figure(1);
                sgtitle('18 _ 10+11, Wndw - 100, Flow -->>')
            end
            if radar == 1
                ax1 = subplot(211);
            else
                ax1 = axes;
            end
            if i == 1
                %surf(ax1,xd,yd,full_target,'EdgeColor','none','FaceColor','interp')
                %hold on
                %plot(ax1,xd, (S.Elevation - bed_A(end:-1:1)),'k',"LineWidth",3)
                xlabel(ax1,'Distance, km')
                ylabel(ax1,'Elevation from sea level, m')
                Z1 = full_target;
                Elevation = S.Elevation;
                bed_fix_A = Elevation - bed_A(end:-1:1);
            end
            if i == 2
                cut = find(~isnan(full_target(:,1)));
                cut = cut(1);
                Z1 = [nan(cut,size(Z1,2)); Z1];
                Z2 = vertcat(full_target,nan(cut,size(full_target,2)));
                Elevation = [S.Elevation Elevation];
                bed_fix_B = S.Elevation - bed_B(end:-1:1);
                %surf(ax1,xd+shift_amount,yd,full_target,'EdgeColor','none','FaceColor','interp')
                %view(180,90)
                %imagesc(ax1,xd,yd,full_target)
            end
            
            %plot(ax1,xd,S.Elevation,'k','LineWidth',3)
            %hold on
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
        Data_size = size(Z1);
        ydd = linspace(max(yd),min(yd),Data_size(1));
        surf(xd(2:end),ydd,[Z2 Z1],'FaceColor','interp','EdgeColor','none')
        view(180,90)
        shading interp
        clim([-1,1])
        colormap(ax1,flipud(CustomColormap))
        cbar = colorbar(ax1);
        set(ax1, 'YDir', 'reverse');
        cbar.Label.String = 'Reflector Slope, Degrees';
        hold on
        plot(xd(2:end),Elevation,'k','LineWidth',3)
        plot(xd(2:end), [bed_fix_B bed_fix_A] ,'k',"LineWidth",3)
    else
        for i = [A B C]
        end
    end
end