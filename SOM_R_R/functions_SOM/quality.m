function q = quality(nn, nn2, x,y,cx,cy)
%QUALITY calculates the "quality" of the retinotopic map
% Description
% - This function calculates the "quality" of the retinotopic map produced
%     - arrays are created that hold the 'x' and 'y' coordinates for a
%     perfect retinotopic map. These are termed 'perfect_x' and 'perfect_y'
%     respectively
%     - The pythagorean distance between the ideal centre of mass and the
%     actual centre of mass of each tectal neuron is calculated and noramlised
%     to the maximum possible distance (sqrt(XT**2+YT**2))
%     - These are summed and divided by the number of tectal neurons
%     - 1 minus this number gives the quality of the map (higher value ->
%     higher quality)
%     
% INPUT
%     XT - integer x dimension of the tectal sheet (number of neurons)
%     YT - integer y dimension of the tectal sheet (number of neurons)
%     XR - integer x dimension of the retinal sheet (number of neurons)
%     YR - integer y dimension of the retinal sheet (number of neurons)
%     COM_X - array of the 'x' coordinates of the tectal neurons' centres
%     of mass
%     COM_Y - array of the 'y' coordinates of the tectal neurons' centres
%     of mass
%
% OUTPUT
%
%     q - integer representing the quality of the retinotopic map produced
%
% COMMENTS  
%
% Adapted from the Python version of George Christopher Baxter at
% https://github.com/geobax/correlated_activity_76
%
% VERSION: Last 2018-05-18  First 2018-05-18
%
% AUTHOR: SSP
% ------------------------------------------------------------------------

x = (x - 0.5)/nn; % The center of each square
y = (y - 0.5)/nn;

maxx = max(x(:));
maxy = max(y(:));
minx = min(x(:));
miny = min(y(:));
maxd = sqrt((maxx-minx)^2+ (maxy-miny)^2);

d2 = sqrt((x'-cx).^2 + (y'-cy).^2);

d2n = d2/maxd;

q = sum(d2n(:))/nn2;

q = 1 - q;

end

