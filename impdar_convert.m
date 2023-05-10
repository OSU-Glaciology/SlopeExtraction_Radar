    function [file_name] = impdar_convert(impdar_mat_path)
%IMPDAR_CONVERT Summary of this function goes here
% Converts data format outputted by ImpDAR to pseudo-CReSIS, just enough
% for SlopeAnalisys
% by Kirill Ivanov

if isstr(impdar_mat_path) == 1
    S = load(impdar_mat_path);
    Time = S.travel_time'/1e+6; %transposed travel_time and converted to seconds (from us) 
    Latitude = S.lat;
    Longitude = -1*S.long; % *1 has been fixed in impdar fork in OSUG, but needs to be here for right now
    Elevation = S.elev;
    Data = single(S.data);
    data_x = S.dist/1000;  %% S assignment is needed for this line, as there is dist function that gets called first
    Data(isnan(Data)) = 0;
    Surface = ones(size(data_x));
    if all(isnan(S.picks.samp2(2,:))) == 1
        bed = S.picks.samp2(1,:);
    else
        bed = S.picks.samp2(2,:); 
    end
    bed = bed*S.picks.pickparams.dt;
    %naming a new restructured .mat file
    %match = wildcardPattern + "/"; isn't in old matlab
    %file_wo_path = erase(fn,match);
    old = '.mat';
    new = '_proc_cresis.mat';
    %file_name = replace(S.fn,old,new);
    file_name = 'proc_cresis.mat';
    save('proc_cresis.mat','Time',"bed","Data","Elevation","Longitude","Latitude","data_x","Surface")
end

