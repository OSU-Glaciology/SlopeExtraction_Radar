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
    if filter == 1
        filt_data = imgaussfilt(full_target,3);
        surf(xd,yd,filt_data,'EdgeColor','none','FaceColor','interp')
    else
        surf(xd,yd,full_target,'EdgeColor','none','FaceColor','interp')
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
                Elevation = S.Elevation;
                D1 = S.Data;
                bed_fix_A = Elevation - bed_A(end:-1:1);
            elseif i == 2
                cut = find(~isnan(S.Data(:,end)));
                cut = cut(1);
                Z1 = [nan(cut,size(Z1,2)); Z1];
                D1 = [nan(cut,size(D1,2)); D1];
                Z2 = vertcat(full_target,nan(cut,size(full_target,2)));
                D2 = vertcat(S.Data,nan(cut,size(S.Data,2)));
                bed_fix_B = S.Elevation - bed_B(end:-1:1);
            end
        end
        Data_size = size(Z1);
        ydd = linspace(max(yd),min(yd),Data_size(1));
        if filter == 1
            filt_data = imgaussfilt([Z2 Z1],3);
            surf(xd(2:end),ydd,filt_data,'FaceColor','interp','EdgeColor','none')
        else
            surf(xd(2:end),ydd,[Z2 Z1],'FaceColor','interp','EdgeColor','none')
        end
        view(180,90)
        shading interp
        clim([-1,1])
        colormap(ax1,CustomColormap)
        cbar = colorbar(ax1);
        set(ax1, 'YDir', 'reverse');
        cbar.Label.String = 'Reflector Slope, Degrees';
        hold on
        plot(xd(2:end),[S.Elevation Elevation],'k','LineWidth',3)
        plot(xd(2:end), [bed_fix_B bed_fix_A] ,'k',"LineWidth",3)
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
                % a - hill on the right, so subset form down, 
                % b subset from the top by the end of a 
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
        if filter ==1
            filt_data = imgaussfilt([Z3 Z2 Z1],3);
            surf(xd(2:end),ydd,filt_data,'FaceColor','interp','EdgeColor','none')
        else
            surf(xd(2:end),ydd,[Z3 Z2 Z1],'FaceColor','interp','EdgeColor','none')
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
        plot(xd(2:end),[E3 E2 E1],'k','LineWidth',3)
        plot(xd(2:end), [bed_fix_C bed_fix_B bed_fix_A] ,'k',"LineWidth",3)
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