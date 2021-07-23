new_hsv = rgb2hsv(frame);
new = new_hsv(:,:,1);

% Removes noise from hue data
new(new > 0.75) = 0;
new(new < 0) = 0;

[height,width] = size(new);
delta = zeros(size(new));

for x = 1:width
    for y = 1:height
         diff = new(y,x);

         % Uses the calibration curve from 
         % CompressionForceCalibration, but will need to recalculate model
         % for different material setups and strain vs. pressure
         delta(y,x) = feval(fittedmodel,diff);
    end
end

%Plot surface mesh
[height_smoothed,width_smoothed] = size(delta);
[X,Y] = meshgrid(1:width_smoothed,1:height_smoothed);
s = surf(X,Y,delta);
colormap(cmap);
s.EdgeColor = 'none';