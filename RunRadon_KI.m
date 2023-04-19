%Code to run RollingRadon (Nick Holschuh)
%Georgia Carroll 
%3/8/2022

%RollingRadon says: 

% The inputs are as follows:
%
% data_x_or_filename - Either the X-axis or the name of a CReSIS flight
%                      file (x-axis in distance)
% data_y - values for the Y-axis ***(ignored if filename provided)
%          this can be either a twtt or depth
% Data - the data raster         ***(ignored if filename provided)
% window - this defines the size of the rolling window
% angle_thresh - this value is the maximum slope useable;
% plotter - 1, generates the debug plots
% [surface_bottom] - a vector containing the surface and bottom picks
% [movie_flag] - 1, Records the debug plots (must have plotter == 1)
% [max_frequency] - this sets the scale for interpolation, based on the
%                 highest f	equency of interest in the data. Can induce
%                 memory problems, and not required.

clear all
close all

%KI RunRadon for Europa with data processed with ImpDAR. 
window = 201;
angle_thresh = [6 6];
plotter = 0;
movie_flag = 0;
load(impdar_convert('19_11_hill.mat'))
%load('proc_cresis.mat')
%surface_bottom=[surface_elev(end)+1e-3-surface_elev;surface_elev(end)+1e-3-bed_elev];
Data = imgaussfilt(real(20*log(Data)));
[r c] = size(Data);
Data = imresize(Data,[round(r/1) round(c/2)],"nearest");
Data(isinf(Data)) = NaN;
[slopegrid_x,slopegrid_y,slopegrid,opt_x,opt_y,opt_angle]=RollingRadon_KI(data_x(1:2:end),Time,Data,window,angle_thresh, ...
    plotter,surface_bottom(1:2:end),movie_flag);



%a=find(new_data(:,end)>0); %bed cutoff
% first cut off bottom of data (Below bed)
%NewData=new_data(a(1):8000,:);
%data_y=surface_elev(end)-depth_axis(a(1):8000); % this needs to be depth, not elevation. 

%addpath('/Users/crosson/Documents/MATLAB/Radargramming/CSARP_standard');
% data_x_or_filename = 'Data_20160512_02_009.mat';
% window = 301;
% angle_thresh = [6 6]; % maybe in degrees? subtracts 5 at one point??
% plotter = 1;
% movie_flag = 0;
% %others
%  load('Data_20160512_02_009.mat');
%  load('BottomQuick0512_02_009.mat');
%[new_data shift_amount new_time] = elevation_shift(Data,Time,Elevation,Surface);
%[new_data shift_amount depth_axis surface_elev bed_elev] = depth_shift(Data,Time,Surface,Elevation,botinter,1);
% [new_data shift_amount depth_axis surface_elev bed_elev] = depth_shift(fliplr(Data),Time,fliplr(Surface),fliplr(Elevation),fliplr(BOTTOM),1);
% a=find(new_data(:,end)>0);




% to get data_x need to convert Lat and Lon into polar
%[x,y]=polarstereo_fwd(Latitude,Longitude,2);
%data_x=distance_vector(x,y,0);



% save 301_66_radon_51602_006.mat

% some plotting functions
                                                                                
% %simple:
% figure; imagesc(slopegrid_x,slopegrid_y,slopegrid)
% xlabel('Distance (meters)')
% ylabel('Depth below highest point (meters)')
% colorbar
% colormap(parula) %you can change this as well
% caxis([-3 3]) % or whatever
% ?
%%
% figure;
[x,y]=meshgrid(slopegrid_x,slopegrid_y);
figure; surf(x,y,-(slopegrid),'EdgeColor','interp');
view(180,90)
colorbar
colormap(redblue)
caxis([-9 9])
hold on
plot(data_x, surface_elev(end)-bed_006+350,'k','LineWidth',3)
% % 
% %% for the radar itself - beware these run kind of slowly
% figure;
%[x,y]=meshgrid(data_x,data_y);
[x,y]=meshgrid(data_x,data_y);
figure; surface(x,y,log10((NewData.^2)),'EdgeColor','interp');
view(180,90)
colormap(bone)
colorbar
caxis([-72 -52]) 
hold on
plot(data_x, surface_elev(end)-bed_006+350,'b','LineWidth',3)

% % % ?
% figure;
% [x,y]=meshgrid(data_x,data_y);
% [x,y]=meshgrid(data_x,data_y);
% figure; f=pcolor(x(1:2:end,1:2:end),y(1:2:end,1:2:end),log10((NewData(1:2:end,1:2:end).^2)))
% shading interp;
% set(gca,'YDir','reverse')
% colormap(bone)
% colorbar
% caxis([5.4 5.75]) 

% %% slices through example
% figure; plot((slopegrid(18:end,620:630)),'r')
% hold on; plot(mean(slopegrid(18:end,620:630)'),'r','LineWidth',5)
% plot((slopegrid(18:end,655:665)),'b')
% hold on; plot(mean(slopegrid(18:end,655:665)'),'b','LineWidth',5)