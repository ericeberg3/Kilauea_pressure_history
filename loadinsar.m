% Open the file in binary mode
filename = 'Data/s20180508_s20180806_corr1.unw.gc';
fid = fopen(filename, 'rb');

% Read the remaining data as 32-bit floats
data = fread(fid, 'float32');
% data = mod(data, 2*pi);
% fclose(fid);
% 
figure; imagesc(reshape(data,4862,4672))
colormap gray;
clim([-10, 10]);

% figure;
% xi = linspace(-8e3, 8e3, 4862);
% yi = linspace(-8e3, 8e3, 4672);
% [X, Y] = meshgrid(xi, yi);
% contourf(X, Y, data, 20, 'LineColor', 'none');
% colormap jet;
colorbar;
% xlabel('Easting (m)', 'FontSize', 14);
% ylabel('Northing (m)', 'FontSize', 14);
% title('InSAR Vertical Displacement', 'FontSize', 16);
% axis equal;