sourceSize = size(test1);
Data_size = size(Data);
targetSize = [Data_size(1)*2-600,Data_size(2)*3];
[X_samples,Y_samples] = meshgrid(linspace(1,sourceSize(2),targetSize(2)), linspace(1,sourceSize(1),targetSize(1)));
source_resized_to_target_size = interp2(test1, X_samples, Y_samples);
figure(1);
ax1 = subplot(211);
imagesc(ax1,source_resized_to_target_size)
colorbar
caxis([-0.4,0.4])
colormap(ax1,redblue(29))
ax2 = subplot(212);
targetSize = [Data_size(1)*2,Data_size(2)*3];
[X_samples,Y_samples] = meshgrid(linspace(1,Data_size(2),targetSize(2)), linspace(1,Data_size(1),targetSize(1)));
D_new = interp2(Data, X_samples, Y_samples);
imagesc(ax2,D_new)
colormap(ax2,'bone')
colorbar
%% combine plot
C = [hill_lake_bump lake_bump];
B = [hill_lake_bump_data lake_bump_data];
figure(1);
ax1 = subplot(211);
imagesc(ax1,C)
colorbar
caxis([-0.4,0.4])
colormap(ax1,redblue(29))
ax2 = subplot(212);
imagesc(ax2,B)
colormap(ax2,'bone')
colorbar