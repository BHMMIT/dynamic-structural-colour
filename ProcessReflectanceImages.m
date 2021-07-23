% Variable parameters (define image folder and number of images)
foldername = "Test set";
imagecount = 1;

% Fixed parameters (sample/lens/screen geometry)
screen_distance_mm = 9;
screen_width_mm = 70;
screen_aperture_dia_mm = 4.76;
lens_na = 0.5;
lens_immersion_n = 1;
lens_distance_mm = 10.6;

% Derived parameters
lens_halfangle_radians = asin(lens_na/lens_immersion_n);
screen_overlap_dia_mm = 2*screen_distance_mm*tan(lens_halfangle_radians);

% Open reference image and draw ROI, clockwise from top left, to get
% moving transformation points
ref = imread(foldername + "/shape_ref.jpg");
imshow(ref);
roi = drawpolygon('StripeColor','y');
moving_points = roi.Position;
close;

% Calculate pixel distances of each side of the trapezoid
ref_roi_top_distance = sqrt((moving_points(2,1) - moving_points(1,1))^2 + (moving_points(2,2) - moving_points(1,2))^2);
ref_roi_bottom_distance = sqrt((moving_points(3,1) - moving_points(4,1))^2 + (moving_points(3,2) - moving_points(4,2))^2);
ref_roi_left_distance = sqrt((moving_points(1,1) - moving_points(4,1))^2 + (moving_points(1,2) - moving_points(4,2))^2);
ref_roi_right_distance = sqrt((moving_points(2,1) - moving_points(3,1))^2 + (moving_points(2,2) - moving_points(3,2))^2);

% Now create manual target transfer points
x_pixels = round((ref_roi_top_distance + ref_roi_bottom_distance)/2);
fixed_points(1,:) = [0;0];
fixed_points(2,:) = [x_pixels;0];
fixed_points(3,:) = [x_pixels;x_pixels];
fixed_points(4,:) = [0;x_pixels];

% Calculate size of each transformed pixel in mm
pixel_width_mm = screen_width_mm/x_pixels;

% Create transform and calculate crop region
transform = fitgeotrans(moving_points,fixed_points,'projective');
[ref_warp,spatial_ref] = imwarp(ref,transform);
fixed_points(:,1) = fixed_points(:,1) - spatial_ref.XWorldLimits(1);
fixed_points(:,2) = fixed_points(:,2) - spatial_ref.YWorldLimits(1);
crop_rectangle = [fixed_points(1,1) fixed_points(1,2) (fixed_points(2,1) - fixed_points(1,1)) (fixed_points(3,2) - fixed_points(1,2))];
ref_warp_crop = imcrop(ref_warp,crop_rectangle);

% Get the center of the screen hole
pinhole_centerpoint = round([x_pixels/2 x_pixels/2]);
pinhole_radius = (screen_aperture_dia_mm/2)/pixel_width_mm;

% Calculate scaling mask
% This accounts for the effect that, the further away a pixel is from the
% center of the screen, the smaller the range of reflection angles that it
% covers.
% NOTE - actually not used to scale images in this version, but worth
% keeping in here
scaling_mask = zeros(size(ref_warp_crop,1),size(ref_warp_crop,2));
for x = 1:size(scaling_mask,2)
    for y = 1:size(scaling_mask,1)
        % Calculate distance to centerpoint
        d = sqrt((x - pinhole_centerpoint(1))^2 + (y - pinhole_centerpoint(2))^2);
        
        % Calculate angle observed by each pixel
        theta = atan((d + 0.5)*pixel_width_mm/screen_distance_mm) - atan((d - 0.5)*pixel_width_mm/screen_distance_mm);
        scaling_mask(y,x) = theta;
    end
end
theta_max = max(max(scaling_mask));
scaling_mask = theta_max./scaling_mask;

% Make pair of circular binary masks and an upper cutoff mask
cutoff_mask = ones(size(ref_warp_crop,1),size(ref_warp_crop,2));
center_mask = ones(size(ref_warp_crop,1),size(ref_warp_crop,2));
overlap_mask = ones(size(ref_warp_crop,1),size(ref_warp_crop,2));
for x = 1:size(scaling_mask,2)
    for y = 1:size(scaling_mask,1)
        % Calculate distance to centerpoint
        d = sqrt((x - pinhole_centerpoint(1))^2 + (y - pinhole_centerpoint(2))^2);
        
        % Set to 0 if inside pinhole
        if (d <= pinhole_radius)
            center_mask(y,x) = 0;
            overlap_mask(y,x) = 0;
        end

        % Set to 0 if outside overlap
        if (d >= (screen_overlap_dia_mm/2)/pixel_width_mm)
            overlap_mask(y,x) = 0;
        end
        
        % Set to 0 if in upper half
        if (y <= pinhole_centerpoint(2))
            cutoff_mask(y,x) = 0;
        end
    end
end

% Get the light reference image and transform it
% NOTE - only needed if background light is not negligible
%light_ref = imwarp(imread(foldername + "/light_ref.jpg"),transform);

% Now loop through all remaining images
for i = 1:imagecount
    % Transform, subtract reference, and crop each image
    current_image = imwarp(imread(foldername + "/" + i + ".jpg"),transform);
    %delta_image_crop = imcrop(current_image - light_ref,crop_rectangle);
    delta_image_crop = imcrop(current_image,crop_rectangle);
    
    % Mask the area of interest
    delta_image_crop_dbl = im2double(delta_image_crop);
    delta_image_crop_dbl(:,:,1) = delta_image_crop_dbl(:,:,1).*cutoff_mask;
    delta_image_crop_dbl(:,:,2) = delta_image_crop_dbl(:,:,2).*cutoff_mask;
    delta_image_crop_dbl(:,:,3) = delta_image_crop_dbl(:,:,3).*cutoff_mask;
    delta_image_crop = im2uint8(delta_image_crop_dbl);
    
    delta_image_crop_dbl_overlap = im2double(delta_image_crop);
    delta_image_crop_dbl_overlap(:,:,1) = delta_image_crop_dbl_overlap(:,:,1).*overlap_mask.*cutoff_mask;
    delta_image_crop_dbl_overlap(:,:,2) = delta_image_crop_dbl_overlap(:,:,2).*overlap_mask.*cutoff_mask;
    delta_image_crop_dbl_overlap(:,:,3) = delta_image_crop_dbl_overlap(:,:,3).*overlap_mask.*cutoff_mask;
    delta_image_crop_overlap = im2uint8(delta_image_crop_dbl_overlap);
    
    % Calculate the total intensity over the whole image
    delta_image_crop_lab = rgb2lab(delta_image_crop);
    delta_image_crop_l = delta_image_crop_lab(:,:,1);
    total_intensity = sum(sum(delta_image_crop_l));
    
    % Calculate the total intensity over the overlap region
    delta_image_crop_overlap_lab = rgb2lab(delta_image_crop_overlap);
    delta_image_crop_overlap_l = delta_image_crop_overlap_lab(:,:,1);
    overlap_intensity = sum(sum(delta_image_crop_overlap_l));
    
    % Add degree overlays to both image
    for theta = 30:20:70
        r = tand(theta)*(screen_distance_mm/pixel_width_mm);
        delta_image_crop = insertShape(delta_image_crop,'circle',[pinhole_centerpoint(1) pinhole_centerpoint(2) r],'LineWidth',3,'Color','white');
    end
    
    % Save files
    imwrite(delta_image_crop,i + ".jpg");
end
