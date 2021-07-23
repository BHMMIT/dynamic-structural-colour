% Needs additional comments and tidying by 7/26/21

%First input forces as an 1xN array called 'forces'

%Import image data
imgfiles = dir('*.tif');
imgcount = numel(imgfiles);
image_data = [];
for i = 1:imgcount
    if i <= 9
        full_image = importdata("ref_0000" + i + ".tif");
    else
        full_image = importdata("ref_000" + i + ".tif");
    end
    % Manually select an ROI for the series of images
    %75_1_8th whole tip
    image_data(:,:,:,i) = full_image(2100:2200,3100:3250,:);
    %757_new whole tip
    %image_data(:,:,:,i) = full_image(2050:2150,3150:3300,:);
    %75_1_8th tiny
    %image_data(:,:,:,i) = full_image(2140:2155,3170:3185,:);
    %75_1_8th pin
    %image_data_pin(:,:,:,i) = full_image(2150,3180,:);
end

%Analyse image regions
k = floor(sqrt(imgcount) + 1);
si_hsv = [];
si_h = [];
N = [];
maxes = [];
averages = [];
edges = [0.04:0.02:1];
contact_area = [];
for i = 1:imgcount
    %Set up different sub-ranges here
    si_hsv(:,:,:,i) = rgb2hsv(image_data(:,:,:,i));
    %For hue average
    %si_hsv(:,:,:,i) = rgb2hsv(image_data(50:70,40:60,:,i));
    temp = si_hsv(:,:,1,i);
    temp(temp>=0.9) = 0;
    temp(temp<=0.04) = 0;
    %temp(si_hsv(:,:,2,i)<0.3) = 0;
    contact_area_pixels(i) = nnz(temp);
    
    %temp_sorted = sort(temp(:),'descend'); 
    %temp_avg = sum(temp_sorted(1:50))/50;
    
    N(i,:) = histcounts(temp,edges);
    %[a,b] = CurveFit(edges(1:end-1),N(i,:));
    %gauss_means(i) = a.b1;
    %maxes(i) = max(max(temp));
    %averages(i) = temp_avg;
    
    figure(2);
    si_h(:,:,i) = temp;
     subplot(k,k,i);
     imshow(si_h(:,:,i));
    
%     figure(2);
%     subplot(k,k,i);
%     imshow(si_hsv(:,:,2,i));
%     
%     figure(3);
%     subplot(k,k,i);
%     imshow(si_hsv(:,:,3,i)./max(max(si_hsv(:,:,3,i))));
end

%Find max contact area
[a,b] = max(contact_area_pixels);
contact_area_pixels(b:end) = a;
contact_area_m2 = contact_area_pixels * 0.00000000019;
%1 pixel = 13.7um x 13.7um = 0.0000137m x 0.0000137m = 0.00000000019 m2
figure(4);
pressure_Pa = (forces-min(forces))./contact_area_m2;
plot(pressure_Pa_lin,hue_average);

%Plot everything
for i = 1:imgcount
    figure(5);
    hold on;
    plot(edges(1:end-1),N(i,:));
    %plot(edges(1:end-1),N(i,:)./max(N));
end

[ver,hor] = size(si_h(:,:,1));
current = [];
hue_curves = [];
%previous = [];
%Plot individual hue runs for each pixel
for i = 1:ver
    for j = 1:hor
        figure(4);
        hold on;
        current(:,1) = si_h(i,j,:);
        flag_0 = 0;
        for k = 5:length(current)
            if current(k,1) == 0
                flag_0 = 1;
            end
        end
        if flag_0 == 0
            %tpc = finddelay(hue_average,current)
            %[c,lags] = xcorr(current,tester);
            %plot(current(tpc+1:end));
            hue_curves(:,end+1) = current(:,1);
            plot(current);
        end
    end
end
hue_average = mean(hue_curves,2);
plot(hue_average);

hue_curves = [];
[ver,hor] = size(si_h(:,:,1));
for i = 1:ver
    for j = 1:hor
        figure(4);
        hold on;
        current(:,1) = si_h(i,j,:);
        hue_curves(:,end+1) = current(:,1);
        plot(current);
    end
end
hue_average = mean(hue_curves,2);
plot(hue_average);

% Nmax = max(N);
% for i = 1:imgcount
%    N(i,:) = N(i,:)./Nmax;
% end
%pcolor(N);