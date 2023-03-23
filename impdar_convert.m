function [file_name] = impdar_convert(impdar_mat_path)
%IMPDAR_CONVERT Summary of this function goes here
% Converts data format outputted by ImpDAR to pseudo-CReSIS, just enough
% for SlopeAnalisys
% by Kirill Ivanov

if isstr(impdar_mat_path) == 1
    load(impdar_mat_path)
    Time = travel_time'/1e+6; %transposed travel_time and converted to seconds (from us) 
    Latitude = lat;
    Longitude = -1*long; % *1 has been fixed in impdar fork in OSUG, but needs to be here for right now
    Elevation = elev;
    Data = data;
    %Data(isnan(Data)) = 0;
    %Surface = picks.samp2(1,:)/1e+6; %Make sure that pick1 is ice surface and pick2 is bottom
    Bottom = picks.samp2(2,:)/1e+6;
    %naming a new restructured .mat file
    match = wildcardPattern + "/";
    file_wo_path = erase(fn,match);
    old = '.mat';
    new = '_proc_cresis.mat';
    file_name = replace(file_wo_path,old,new);
    save(file_name,'Time',"Bottom","Data","Elevation","Longitude","Latitude")
end

