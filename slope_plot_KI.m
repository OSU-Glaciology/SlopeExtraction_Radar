%% Kirill Ivanov 
%run to plot resuslts from RollingRadon slope analysis by Nich Holschun
%% Init
clear all
close all


% plot radargram as a subplot with the slope plot
radar = 1; % 0 - no, 1 - yes 

%choose line
line = 1; % 0 - 19_11, 1 - 18_10-11, 2 - 18_6-7-8

%gaussian filter the slope output
filter = 1; % 0 - no, 1 - yes

%vertical profile of the reflector slope
%returns a separate figure with a profile.
vert_profile = 1; % 0 - no, 1 - yes. 

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
    A = load('18_6_full.mat');
    bed_A = readmatrix("18_6_bot.csv");
    bed_A = bed_A(:,4)';
    B = load("18_7_full.mat");
    bed_B = readmatrix("18_7_bot.csv");
    bed_B = bed_B(:,4)';
    C = load("18_8_full.mat");
    bed_C = readmatrix("18_8_bot.csv");
    bed_C = bed_C(:,4)';

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
    Slope = full_target;
    if filter == 1
        Slope = imgaussfilt(Slope,3);
        surf(xd,yd,Slope,'EdgeColor','none','FaceColor','interp')
    else
        surf(xd,yd,Slope,'EdgeColor','none','FaceColor','interp')
    end
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
    bed_elev = Elevation - bed(end:-1:1);
    plot(ax1,xd, bed_elev,'k',"LineWidth",3)
    area(ax1, xd, bed_elev-10,min(yd),'FaceColor','white', ...
        'EdgeColor','none','FaceAlpha',1)
    ylim(ax1,[min(yd) max(yd)])
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
        figure(1);
        sgtitle('18 _ 10+11, Wndw - 100, Flow <<--')
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
            if radar == 1
                ax1 = subplot(211);
            else
                ax1 = axes;
            end
            if i == 1
                xlabel(ax1,'Distance, km')
                ylabel(ax1,'Elevation from sea level, m')
                Z1 = full_target;
                E1 = S.Elevation;
                D1 = S.Data;
                bed_fix_A = E1 - bed_A(end:-1:1);
            elseif i == 2
                cut = find(~isnan(S.Data(:,end)));
                cut = cut(1);
                Z1 = [nan(cut,size(Z1,2)); Z1];
                D1 = [nan(cut,size(D1,2)); D1];
                Z2 = vertcat(full_target,nan(cut,size(full_target,2)));
                D2 = vertcat(S.Data,nan(cut,size(S.Data,2)));
                E2 = S.Elevation;
                bed_fix_B = E2 - bed_B(end:-1:1);
            end
        end
        Data_size = size(Z1);
        ydd = linspace(max(yd),min(yd),Data_size(1));
        Elevation = [E2 E1];
        bed_elev = [bed_fix_B bed_fix_A];
        Slope = [Z2 Z1];
        if filter == 1
            Slope = imgaussfilt(Slope,3);
            surf(xd(2:end),ydd,Slope,'FaceColor','interp','EdgeColor','none')
        else
            surf(xd(2:end),ydd,Slope,'FaceColor','interp','EdgeColor','none')
        end
        view(180,90)
        shading interp
        clim([-1,1])
        colormap(ax1,CustomColormap)
        cbar = colorbar(ax1);
        set(ax1, 'YDir', 'reverse');
        cbar.Label.String = 'Reflector Slope, Degrees';
        hold on
        plot(xd(2:end),Elevation,'k','LineWidth',3)
        plot(xd(2:end), bed_elev ,'k',"LineWidth",3)
        area(xd(2:end), bed_elev-10,min(ydd),'FaceColor','white', ...
        'EdgeColor','none','FaceAlpha',1)
        ylim(ax1,[min(ydd) max(ydd)])
        % Radargram
        if radar == 1
            ax2 = subplot(212);
            imagesc(ax2,xd(2:end), ydd(end:-1:1), [D2 D1])
            colormap(ax2,'bone')
            colorbar
            set(ax2, 'XDir', 'reverse');
            xlabel(ax2,'Distance, km')
            ylabel(ax2,'Elevation from sea level, m')
            linkaxes([ax1, ax2], 'x','y');
        end
    else
        xd = zeros(1);
        yd = [];
        figure(1);
        sgtitle('18 _ 6+7+8, Wndw - 100, Flow <<--')
        for i = 1:3
            if i == 1
                S = A;
            elseif  i == 2
                S = B;
            else
                S = C;
            end
            %slope
            slp = -S.slopegrid;
            slp(slp == 6) = NaN;
            yd = [yd; max(S.Elevation) - (S.Time/2 * 1.69e8)];
            xd = [xd (S.data_x*1000 + max(xd))];
            %resize
            sourceSize = size(slp);
            if i == 2
                Data_size = size(S.Data) + [2,0];
            else
                Data_size = size(S.Data);
            end
            targetSize = [Data_size(1)-S.window+1,Data_size(2)];
            [X_samples,Y_samples] = meshgrid(linspace(1,sourceSize(2),targetSize(2)), linspace(1,sourceSize(1),targetSize(1)));
            source_resized_to_target_size = interp2(slp, X_samples, Y_samples);
            full_target = [nan((S.window-1)/2,size(source_resized_to_target_size,2)); source_resized_to_target_size];
            full_target = vertcat(full_target,nan((S.window-1)/2,size(source_resized_to_target_size,2)));
            if radar == 1
                ax1 = subplot(211);
            else
                ax1 = axes;
            end
            if i == 1
                xlabel(ax1,'Distance, km')
                ylabel(ax1,'Elevation from sea level, m')
                Z1 = full_target;
                E1 = S.Elevation;
                D1 = S.Data;
                bed_fix_A = E1 - bed_A(end:-1:1);
            elseif i == 2
                Data_size = size(D1);
                temp_size = size(S.Data);
                Z2 = vertcat(full_target,nan(Data_size(1)-temp_size(1),size(full_target,2)));
                D2 = vertcat(S.Data,nan(Data_size(1)-temp_size(1),size(S.Data,2)));
                cut = find(~isnan(D1(:,2)));
                cut = cut(1);
                Z2 = [nan(cut,size(Z2,2)); Z2];
                D2 = [nan(cut,size(D2,2)); D2];
                Z1 = vertcat(Z1,nan(cut+2,size(Z1,2)));
                D1 = vertcat(D1,nan(cut+2,size(D1,2)));
                E2 = S.Elevation;
                bed_fix_B = E2 - bed_B(end:-1:1);
                
            else
                cut = find(~isnan(D2(:,2)));
                cut = cut(1);
                Z3 = [nan(cut,size(full_target,2)); full_target];
                D3 = [nan(cut,size(S.Data,2)); S.Data];
                Data_size = size(D1);
                temp_size = size(D3);
                Z3 = vertcat(Z3,nan(Data_size(1)-temp_size(1),size(Z3,2)));
                D3 = vertcat(D3,nan(Data_size(1)-temp_size(1),size(D3,2)));
                E3 = S.Elevation;
                bed_fix_C = E3 - bed_C(end:-1:1);
            end
        end
        Data_size = size(Z1);
        ydd = linspace(max(yd),min(yd),Data_size(1));
        Elevation = [E3 E2 E1];
        bed_elev = [bed_fix_C bed_fix_B bed_fix_A];
        Slope = [Z3 Z2 Z1];
        if filter ==1
            Slope = imgaussfilt(Slope,3);
            surf(xd(2:end),ydd,Slope,'FaceColor','interp','EdgeColor','none')
        else
            surf(xd(2:end),ydd,Slope,'FaceColor','interp','EdgeColor','none')
        end
        view(180,90)
        shading interp
        clim([-1,1])
        colormap(ax1,CustomColormap)
        cbar = colorbar(ax1);
        set(ax1, 'YDir', 'reverse');
        cbar.Label.String = 'Reflector Slope, Degrees';
        axis tight
        hold on
        plot(xd(2:end),Elevation,'k','LineWidth',3)
        plot(xd(2:end), bed_elev ,'k',"LineWidth",3)
        area(xd(2:end), bed_elev-10,min(ydd),'FaceColor','white', ...
        'EdgeColor','none','FaceAlpha',1)
        ylim(ax1,[min(ydd) max(ydd)])
        % Radargram
        if radar == 1
            ax2 = subplot(212);
            imagesc(ax2,xd(2:end), ydd(end:-1:1), [D3(1:1:end-2,:) D2 D1(1:1:end-2,:)])
            colormap(ax2,'bone')
            colorbar
            set(ax2, 'XDir', 'reverse');
            xlabel(ax2,'Distance, km')
            ylabel(ax2,'Elevation from sea level, m')
            linkaxes([ax1, ax2], 'x','y');
        end
    end
end
%% vertical Profile 
if vert_profile == 1
    if filter == 0
        disp('Recommend to look at vertical slope profile with filtering.')
    end
    prompt = 'Click on the slope graph for x input: ';
    title = 'vertical slope profile';
    xval = inputdlg(prompt,title);
    xval = str2double(xval{1});
    if line ~= 0
        xd = xd(2:end);
        yd = ydd;
    end
    [~, idx] = min(abs(xd - xval));
    hold on;
    plot(xd(idx), yd(idx), 'ro');
    [xclick, ~] = ginput(1);
    disp(['X value for vertical profile: ', num2str(xclick),'km.']);
    hold on 
    plot(ax1,[xclick xclick],[max(bed_elev) min(Elevation)],'Color','magenta','LineWidth',2)
    hold off
    figure(2)
    sgtitle(['Verical Profile of Reflector Slope at: ', num2str(xclick),'km.'])
    ax3 = axes;
    xval = find((xd - xclick)>0);
    xval = xval(1);
    yval = find((ydd - max(bed_elev))>0);
    x = Slope(1:max(yval),xval);
    y = ydd(1:max(yval));
    idx = (x < 0);
    plot(x(idx),y(idx),'r*',x(~idx),y(~idx),'b*')
    hold on
    k = isfinite(x);
    trend = polyfit(x(k), y(k), 2);
    px = [min(x) max(x)];
    py = polyval(trend, px);
    plot(ax3,px, py, 'LineWidth', 2,'LineStyle',':','Color','black');
    xlabel(ax3,'Reflector Slope, deg')
    ylabel(ax3,'Elevation from sea level, m') 
end