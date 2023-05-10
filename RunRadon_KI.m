% Code to run RollingRadon (Nick Holschuh)
% Kirill Ivanov
% includes graphing
clear global
close all
%% Set-up
% General
    %choose line
    line = 2; % 0 - 19_11, 1 - 18_10-11, 2 - 18_6-7-8

% Radon setup
    SKIP_RADON = 0; %skip radon process for plotting only with existing mat files
    window = 61; % pixel window 
    angle_thresh = [5 5]; % reflector slope limits in degrees
    plotter = 0; % 1 - plotter for slope rolling
    movie_flag = 0;

%Plotting setup
    % plot radargram as a subplot with the slope plot
    radar = 1; % 0 - no, 1 - yes 
    %gaussian filter the slope output
    filter = 0; % 0 - no, 1 - yes
    %vertical profile of the reflector slope
    %returns a separate figure with a profile.
    vert_profile = 1; % 0 - no, 1 - yes. 
        wind_vert_prof = 250; % x window for vertical profile to average
%% Load
load('CustomColormap.mat')
if SKIP_RADON == 0
    if line == 0
        load("19_11_full.mat")
        layer = readmatrix('19_11.csv');
        layer = layer(:,(4:11))';
    elseif line == 1
        A = load(impdar_convert('18_10_flat.mat'));
        %layer_A = readmatrix("18_10.csv");
        %layer_A = layer_A(:,(4:10))';
        B = load(impdar_convert('18_11_flat.mat'));
        %layer_B = readmatrix("18_11.csv");
        %layer_B = layer_B(:,(4:10))';
        Time = [A.Time; B.Time];
        Latitude = [A.Latitude B.Latitude];
        Longitude = [A.Longitude B.Longitude];
        Elevation = [A.Elevation B.Elevation];
        Data = [A.Data B.Data];
        Surface = [A.Surface B.Surface];
        bed = [A.bed B.bed];
        ydd = min(Time):8.00e-9:max(Time);
        Time = ydd';
    else
        C = load(impdar_convert('18_6_flat.mat'));
        % layer_A = readmatrix("18_6.csv");
        % layer_A = layer_A(:,(4:8))';
        B = load(impdar_convert('18_7_flat.mat'));
        % layer_B = readmatrix("18_7.csv");
        % layer_B = layer_B(:,(4:8))';
        A = load(impdar_convert('18_8_flat.mat'));
        % layer_C = readmatrix("18_8.csv");
        % layer_C = layer_C(:,(4:8))';
        Time = [A.Time; B.Time; C.Time];
        Latitude = [A.Latitude B.Latitude C.Latitude];
        Longitude = [A.Longitude B.Longitude C.Longitude];
        Elevation = [A.Elevation B.Elevation C.Elevation];
        Data = [A.Data B.Data C.Data];
        Surface = [A.Surface B.Surface C.Surface];
        bed = [A.bed B.bed C.bed];
        ydd = min(Time):8.00e-9:max(Time);
        Time = ydd';
    end
else
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
end
%% Run Radon 
if SKIP_RADON == 0
    [new_data,shift_amount,depth_axis,surface_elev,bed_elev] = depth_shift(Data,Time,Surface,Elevation,bed,1);
    Data = imgaussfilt(real(20*log(new_data)));
    data_y=surface_elev(end)-depth_axis; 
    [x,y]=polarstereo_fwd(Latitude,Longitude,2);
    data_x=distance_vector(x,y,0);
    surface_bottom=[surface_elev(end)+1e-3-surface_elev;surface_elev(end)+1e-3-bed_elev];
    [slopegrid_x,slopegrid_y,slopegrid,opt_x,opt_y,opt_angle]=RollingRadon(data_x,data_y,Data,window,angle_thresh, ...
        plotter,surface_bottom,movie_flag);
end
%% Graphing
slp = -slopegrid;
slp(slp == angle_thresh(1)) = NaN;
yd = max(surface_elev) - Time/2 * 1.69e8;
xd = data_x;
sourceSize = size(slp);
Data_size = size(Data);
targetSize = [Data_size(1)-window+1,Data_size(2)];
[X_samples,Y_samples] = meshgrid(linspace(1,sourceSize(2),targetSize(2)), linspace(1,sourceSize(1),targetSize(1)));
source_resized_to_target_size = interp2(slp, X_samples, Y_samples);
full_target = [nan((window-1)/2,size(source_resized_to_target_size,2)); source_resized_to_target_size];
full_target = vertcat(full_target,nan((window-1)/2,size(source_resized_to_target_size,2)));

figure(1);
if line == 0
    sgtitle('19 _ 11, Wndw - 100, Flow -->>')
elseif line == 1
    sgtitle('18 _ 10+11, Wndw - 100, Flow <<--')
else
    sgtitle('18 _ 6+7+8, Wndw - 100, Flow -->>')
end
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
clim([-angle_thresh(1),angle_thresh(1)])
colormap(ax1,flipud(CustomColormap))
cbar = colorbar(ax1);
set(ax1, 'YDir', 'reverse');
cbar.Label.String = 'Reflector Slope, Degrees';
xlabel(ax1,'Distance, km')
ylabel(ax1,'Elevation from sea level, m')
hold on
plot(ax1,xd,surface_elev,'k','LineWidth',3)
%layer_elev = Elevation - layer(:,end:-1:1);
plot(ax1,xd, bed_elev,'k',"LineWidth",3)
area(ax1, xd, bed_elev-10,min(yd),'FaceColor','white', ...
    'EdgeColor','none','FaceAlpha',1)
% for i = 2:7
%     plot(ax1,xd,layer_elev(i,:),"Color",[0.7 0.7 0.7],'LineWidth',1.5)
% end
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

% clearvars -except ax1 ax2 xd yd ydd D1 D2 D3 Slope Elevation ...
%     layer_elev vert_profile line filter cbar CustomColormap wind_vert_prof
%% vertical Profile 
if vert_profile == 1
    if filter == 0
        disp('Recommend to look at vertical slope profile with filtering.')
    end
    prompt = 'Click on the slope graph for x input: ';
    title = 'vertical slope profile';
    xval = inputdlg(prompt,title);
    xval = str2double(xval{1});
    [~, idx] = min(abs(xd - xval));
    hold on;
    plot(xd(idx), yd(idx), 'ro');
    [xclick, ~] = ginput(1);
    disp(['X value for vertical profile: ', num2str(xclick)/1000,'km.']);
    hold on 
    plot(ax1,[xclick xclick],[max(bed_elev) min(surface_elev)],'Color','magenta','LineWidth',2)
    hold off
    figure(2)
    sgtitle(['Verical Profile of Reflector Slope at: ', num2str(xclick)/1000,'km.'])
    ax3 = axes;
    xval = find((xd - xclick)>0);
    xval = xval(1);
    yval = find((yd - max(bed_elev))>0);
    x = Slope(1:max(yval),xval-wind_vert_prof:xval);
    y = yd(1:max(yval));
    x_avg = mean(x,2);
    idx = (x_avg < 0);
    plot(ax3,x_avg(idx),y(idx),'r*',x_avg(~idx),y(~idx),'b*')
    %plot(x_avg,y)
    %colormap(ax3,CustomColormap)
    clim([-angle_thresh(1),angle_thresh(1)])
    hold on
    k = ~isnan(x_avg);
    p = polyfit(x_avg(k), y(k), 1);
    f = polyval(p,x_avg);
    plot(ax3,x_avg, f, 'LineWidth', 2,'LineStyle',':','Color','black');
    xlabel(ax3,'Reflector Slope, deg')
    ylabel(ax3,'Elevation from sea level, m') 
    xlim(ax3,[-angle_thresh(1),angle_thresh(1)])
    ylim(ax3,[min(yd) max(yd)])
end

fig = figure(4);
ax1 = subplot(1,3,1,'Parent',fig);
clear title
Slope_f = imgaussfilt(Slope,2);
surf(ax1,xd,yd,Slope,'FaceColor','interp','EdgeColor','none')
if line == 1
    title('I9-C9 Line Data Slope Extraction','FontSize',16)
elseif line == 2
    title('A9-012 Line Data Slope Extraction','FontSize',16)
else
    title('I10-C10 Line Data Slope Extration','FontSize',16)
end
view(180,90)
shading interp
clim([-angle_thresh(1),angle_thresh(1)])
colormap(ax1,CustomColormap)
cbar = colorbar(ax1,"eastoutside");
set(ax1, 'YDir', 'reverse');
cbar.Label.String = 'Reflector Slope, Degrees';
set(cbar,'Position',cbar.Position + [0.4 0.15 0 -0.3])
xlabel(ax1,'Distance, km')
hold on
ylim(ax1,[min(yd) max(yd)])
box on
% a_x = [xclick xclick]; 
% a_y = [0.85 0.85];      
% annotation('line',a_x,a_y,'String',' Flow Direction ','FontSize',9,'Linewidth',2)
if line == 1
    a_x = [0.5 0.45]; 
    a_y = [0.86 0.86];      
    annotation('textarrow',a_x,a_y,'String',' Flow Direction ','FontSize',11,'Linewidth',3)
elseif  line == 2
    a_x = [0.45 0.5]; 
    a_y = [0.83 0.83];      
    annotation('textarrow',a_x,a_y,'String',' Flow Direction ','FontSize',11,'Linewidth',3)   
end
hold off

%radar
ax = subplot(1,3,2);
imagesc(ax,xd,yd,Data,'AlphaData',0.5)
hold on
h(1) = plot(ax,[xclick xclick],[max(bed_elev) min(surface_elev)],'Color','black', ...
    'LineWidth',1);
h(2) = plot(ax,[xclick-(wind_vert_prof*0.003) xclick-(wind_vert_prof*0.003)], ...
    [max(bed_elev) min(surface_elev)],'Color','black', ...
    'LineWidth',1);
h(3) = plot(ax,[xclick xclick-(wind_vert_prof*0.003)],[max(bed_elev) max(bed_elev)], ...
    'Color','black','LineWidth',1);
h(3) = plot(ax,[xclick xclick-(wind_vert_prof*0.003)],[min(surface_elev) min(surface_elev)], ...
    'Color','black','LineWidth',1);
uistack(h,'top')
plot(xd,surface_elev,'k','LineWidth',3)
plot(xd, bed_elev,'k',"LineWidth",3)
area(xd, bed_elev-10,min(yd),'FaceColor','white', ...
'EdgeColor','none','FaceAlpha',1)
% max_i = size(layer_elev(:,1));
% for i = 2:max_i(1)
%     plot(xd,layer_elev(i,:),'Color',[.7 .7 .7],'LineWidth',1.5)
% end
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
xlim(ax2,[-angle_thresh(1),angle_thresh(1)])
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