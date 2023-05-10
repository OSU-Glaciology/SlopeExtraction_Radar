%% Kirill Ivanov 
%run to plot resuslts from RollingRadon slope analysis by Nich Holschun
%% Init
clear global
%close all


% plot radargram as a subplot with the slope plot
radar = 1; % 0 - no, 1 - yes 

%choose line
line = 1; % 0 - 19_11, 1 - 18_10-11, 2 - 18_6-7-8

%gaussian filter the slope output
filter = 0; % 0 - no, 1 - yes

%vertical profile of the reflector slope
%returns a separate figure with a profile.
vert_profile = 1; % 0 - no, 1 - yes. 
    wind_vert_prof = 150; % x window for vertical profile to average

%% Load
load('CustomColormap.mat')
if line == 0
    load("19_11_full.mat")
    layer = readmatrix('19_11.csv');
    layer = layer(:,(4:11))';
elseif line == 1
    A = load("18_10_full_2.mat");
    layer_A = readmatrix("18_10.csv");
    layer_A = layer_A(:,(4:10))';
    B = load("18_11_full_2.mat");
    layer_B = readmatrix("18_11.csv");
    layer_B = layer_B(:,(4:10))';
else
    A = load('18_6_full.mat');
    layer_A = readmatrix("18_6.csv");
    layer_A = layer_A(:,(4:8))';
    B = load("18_7_full.mat");
    layer_B = readmatrix("18_7.csv");
    layer_B = layer_B(:,(4:8))';
    C = load("18_8_full.mat");
    layer_C = readmatrix("18_8.csv");
    layer_C = layer_C(:,(4:8))';
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
    layer_elev = Elevation - layer(:,end:-1:1);
    plot(ax1,xd, layer_elev(1,:),'k',"LineWidth",3)
    area(ax1, xd, layer_elev(1,:)-10,min(yd),'FaceColor','white', ...
        'EdgeColor','none','FaceAlpha',1)
    for i = 2:7
        plot(ax1,xd,layer_elev(i,:),"Color",[0.7 0.7 0.7],'LineWidth',1.5)
    end
    ylim(ax1,[min(yd) max(yd)])
    % Radargram
    if radar == 1
        ax2 = subplot(212);
        imagesc(ax2,xd, yd, Data)
        colormap(ax2,'bone')
        set(ax2, 'XDir', 'reverse');
        set(ax2,'YDir','normal');
        xlabel(ax2,'Distance, km')
        ylabel(ax2,'Elevation from sea level, m')
        pos1 = get(ax1,'Position');
        pos2 = get(ax2,'Position');
        pos1(3) = pos2(3);
        set(ax1,'Position',pos1)
        axis tight
        linkaxes([ax1, ax2], 'xy');
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
                layer_A = E1 - layer_A(:,end:-1:1);
            elseif i == 2
                cut = find(~isnan(S.Data(:,end)));
                cut = cut(1);
                Z1 = [nan(cut,size(Z1,2)); Z1];
                D1 = [nan(cut,size(D1,2)); D1];
                Z2 = vertcat(full_target,nan(cut,size(full_target,2)));
                D2 = vertcat(S.Data,nan(cut,size(S.Data,2)));
                E2 = S.Elevation;
                layer_B = E2 - layer_B(:,end:-1:1);
            end
        end
        Data_size = size(Z1);
        ydd = linspace(max(yd),min(yd),Data_size(1));
        Elevation = [E2 E1];
        layer_elev = [layer_B layer_A];
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
        plot(xd(2:end), layer_elev(1,:) ,'k',"LineWidth",3)
        area(xd(2:end), layer_elev(1,:)-10,min(ydd),'FaceColor','white', ...
        'EdgeColor','none','FaceAlpha',1)
        for i = 2:6
            plot(xd(2:end),layer_elev(i,:),'Color',[.7 .7 .7],'LineWidth',1.5)
        end
        ylim(ax1,[min(ydd) max(ydd)])
        % Radargram
        if radar == 1
            ax2 = subplot(212);
            imagesc(ax2,xd(2:end), ydd, [D2 D1])
            colormap(ax2,'bone')
            set(ax2,'YDir','normal');
            set(ax2, 'XDir', 'reverse');
            xlabel(ax2,'Distance, km')
            ylabel(ax2,'Elevation from sea level, m')
            pos1 = get(ax1,'Position');
            pos2 = get(ax2,'Position');
            pos1(3) = pos2(3);
            set(ax1,'Position',pos1)
            axis tight
            linkaxes([ax1, ax2], 'xy');
        end
    else
        xd = zeros(1);
        yd = [];
        figure(1);
        sgtitle('18 _ 6+7+8, Wndw - 100, Flow -->>')
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
                layer_A = E1 - layer_A(:,end:-1:1);
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
                layer_B = E2 - layer_B(:,end:-1:1);
                
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
                layer_C = E3 - layer_C(:,end:-1:1);
            end
        end
        Data_size = size(Z1);
        ydd = linspace(max(yd),min(yd),Data_size(1));
        Elevation = [E3 E2 E1];
        layer_elev = [layer_C layer_B layer_A];
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
        plot(xd(2:end), layer_elev(1,:) ,'k',"LineWidth",3)
        area(xd(2:end), layer_elev(1,:)-10,min(ydd),'FaceColor','white', ...
        'EdgeColor','none','FaceAlpha',1)
        for i = 2:5
            plot(xd(2:end),layer_elev(i,:),'Color',[.7 .7 .7],'LineWidth',1.5)
        end
        ylim(ax1,[min(ydd) max(ydd)])
        % Radargram
        if radar == 1
            ax2 = subplot(212);
            imagesc(ax2,xd(2:end), ydd, [D3(1:1:end-2,:) D2 D1(1:1:end-2,:)])
            colormap(ax2,'bone')
            set(ax2,'YDir','normal');
            set(ax2, 'XDir', 'reverse');
            xlabel(ax2,'Distance, km')
            ylabel(ax2,'Elevation from sea level, m')
            pos1 = get(ax1,'Position');
            pos2 = get(ax2,'Position');
            pos1(3) = pos2(3);
            set(ax1,'Position',pos1)
            axis tight
            linkaxes([ax1, ax2], 'xy');
        end
    end
end

clearvars -except ax1 ax2 xd yd ydd D1 D2 D3 Slope Elevation ...
    layer_elev vert_profile line filter cbar CustomColormap wind_vert_prof
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
    plot(ax1,[xclick xclick],[max(layer_elev(1,:)) min(Elevation)],'Color','magenta','LineWidth',2)
    hold off
    figure(2)
    sgtitle(['Verical Profile of Reflector Slope at: ', num2str(xclick),'km.'])
    ax3 = axes;
    xval = find((xd - xclick)>0);
    xval = xval(1);
    yval = find((yd - max(layer_elev(1,:)))>0);
    x = Slope(1:max(yval),xval-wind_vert_prof:xval);
    y = yd(1:max(yval));
    x_avg = mean(x,2);
    idx = (x_avg < 0);
    plot(ax3,x_avg(idx),y(idx),'r*',x_avg(~idx),y(~idx),'b*')
    %plot(x_avg,y)
    %colormap(ax3,CustomColormap)
    clim([-0.5 0.5])
    hold on
    k = ~isnan(x_avg);
    p = polyfit(x_avg(k), y(k), 1);
    f = polyval(p,x_avg);
    plot(ax3,x_avg, f, 'LineWidth', 2,'LineStyle',':','Color','black');
    xlabel(ax3,'Reflector Slope, deg')
    ylabel(ax3,'Elevation from sea level, m') 
    xlim(ax3,[-1 1])
    ylim(ax3,[min(yd) max(yd)])
end
%% figure plotting
fig = figure(4);
ax1 = subplot(1,3,1,'Parent',fig);
clear title
Slope_f = imgaussfilt(Slope,2);
surf(ax1,xd,ydd,Slope,'FaceColor','interp','EdgeColor','none')
title('I9-C9 Line Data Slope Extraction','FontSize',16)
view(180,90)
shading interp
clim([-1,1])
colormap(ax1,CustomColormap)
cbar = colorbar(ax1,"eastoutside");
set(ax1, 'YDir', 'reverse');
cbar.Label.String = 'Reflector Slope, Degrees';
set(cbar,'Position',cbar.Position + [0.4 0.15 0 -0.3])
xlabel(ax1,'Distance, km')
hold on
ylim(ax1,[min(ydd) max(ydd)])
box on
% a_x = [xclick xclick]; 
% a_y = [0.85 0.85];      
% annotation('line',a_x,a_y,'String',' Flow Direction ','FontSize',9,'Linewidth',2)
a_x = [0.5 0.45]; 
a_y = [0.86 0.86];      
annotation('textarrow',a_x,a_y,'String',' Flow Direction ','FontSize',11,'Linewidth',3)
hold off

%radar
ax = subplot(1,3,2);
if line == 1
    imagesc(ax,xd,ydd,[D2 D1],'AlphaData',0.5)
elseif line == 2
    imagesc(ax,xd,ydd,[D3(1:1:end-2,:) D2 D1(1:1:end-2,:)],'AlphaData',0.5)
else
    imgesc(ax,xd,ydd,D1,'AlphaData',0.5)
end
hold on
h(1) = plot(ax,[xclick xclick],[max(layer_elev(1,:)) min(Elevation)],'Color','black', ...
    'LineWidth',1);
h(2) = plot(ax,[xclick-(wind_vert_prof*0.003) xclick-(wind_vert_prof*0.003)], ...
    [max(layer_elev(1,:)) min(Elevation)],'Color','black', ...
    'LineWidth',1);
h(3) = plot(ax,[xclick xclick-(wind_vert_prof*0.003)],[max(layer_elev(1,:)) max(layer_elev(1,:))], ...
    'Color','black','LineWidth',1);
h(3) = plot(ax,[xclick xclick-(wind_vert_prof*0.003)],[min(Elevation) min(Elevation)], ...
    'Color','black','LineWidth',1);
uistack(h,'top')
plot(xd,Elevation,'k','LineWidth',3)
plot(xd, layer_elev(1,:) ,'k',"LineWidth",3)
area(xd, layer_elev(1,:)-10,min(ydd),'FaceColor','white', ...
'EdgeColor','none','FaceAlpha',1)
max_i = size(layer_elev(:,1));
for i = 2:max_i(1)
    plot(xd,layer_elev(i,:),'Color',[.7 .7 .7],'LineWidth',1.5)
end
linkaxes([ax1 ax],'xy')
ax.Visible = 'off';
ax.XTick = [];
ax.YTick = [];
colormap(ax,'bone');
set(ax,'YDir','normal');
set(ax, 'XDir', 'reverse');
pos1 = get(ax1,'Position');
pos1 = pos1 + [0.01 0 0.3 0];
set([ax1 ax],'Position',pos1);
hold off

%trend
ax2 = subplot(1,3,3);
plot(ax2,x_avg,y)
title('Verical Profile of Reflector Slope','FontSize',16)
xlabel(ax2,'Reflector Slope, deg')
xlim(ax2,[-1 1])
ylim(ax2,[min(yd) max(yd)])
set(gca,'yticklabel',[])
box on
h = axes(fig,'visible','off');
pos2 = get(ax2,'Position');
pos2 = pos2 + [0.05 0 -0.01 0];
set(ax2,'Position',pos2)
axis tight
h.YLabel.Visible = 'on';
ylabel(h,'Elevation from sea level, m','FontWeight','bold');
linkaxes([ax1, ax2], 'y');