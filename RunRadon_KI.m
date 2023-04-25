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
window = 101;
angle_thresh = [6 6];
plotter = 0;
movie_flag = 0;
load(impdar_convert('18_11_bot_cut.mat'))
%surface_bottom=[surface_elev(end)+1e-3-surface_elev;surface_elev(end)+1e-3-bed_elev];
Data = imgaussfilt(real(20*log(Data)));
%[r c] = size(Data);
%Data = imresize(Data,[round(r/1) round(c/2)],"nearest");
Data(isinf(Data)) = NaN;
[slopegrid_x,slopegrid_y,slopegrid,opt_x,opt_y,opt_angle]=RollingRadon_KI(data_x,Time,Data,window,angle_thresh, ...
    plotter,surface_bottom,movie_flag);

save('18_11_full.mat')