% This file plots measured spectra, taken from a sample being stretched,
% as well as converting those spectra into points in the CIE 1931 color
% space

% The inputs to this are:
% WavelengthData - the array of wavelengths the spectrometer measures at
% (Ocean Optics Maya 2000)
% SpectrumData - the actual measurements at those wavelengths
% ForceData - array of measured force data at each point of strain
% StartLengthText - initial length of sample, as a string

% Number of datapoints
num_datapoints = size(WavelengthData,2);

% Spectral range of interest, defined by user
wavelength_lower = 450;
wavelength_upper = 650;

% Find the closest values in the spectrometer wavelength data
[~,wavelength_index_lower] = min(abs(WavelengthData(:,1)-wavelength_lower));
[~,wavelength_index_upper] = min(abs(WavelengthData(:,1)-wavelength_upper));
wavelength_array = WavelengthData(wavelength_index_lower:wavelength_index_upper,1);

% Prepare color matching functions using https://www.mathworks.com/matlabcentral/fileexchange/7021-spectral-and-xyz-color-functions
[cie_l,cie_x,cie_y,cie_z] = colorMatchFcn('1931_full');

% Setup interpolated color matching functions
interpolated_cie_x = interp1(cie_l,cie_x,wavelength_array);
interpolated_cie_y = interp1(cie_l,cie_y,wavelength_array);
interpolated_cie_z = interp1(cie_l,cie_z,wavelength_array);

% Convert force data from string to double
ForceData = str2double(ForceData);

% Convert start length from string to double
StartLength = str2double(StartLengthText);

% Initialize CIE variables
coordinate_x = 0;
coordinate_y = 0;
coordinate_rgb = [0 0 0];

% Initialize spectra plot
figure(1);
hold on

for i = 1:num_datapoints
    % Get spectrum for current datapoint
    spectral_array = SpectrumData(wavelength_index_lower:wavelength_index_upper,i);
    
    % Calculate strain value
    current_strain = abs(PositionData(1,i) - PositionData(1,1))/StartLength
    strain_array = current_strain*ones(size(wavelength_array,1),1);
    
    % Integrate reflectivity
    trapz(spectral_array)
    
    % Calculate CIE xy coordinates to correctly color each individual curve
    integral_X = sum(spectral_array.*interpolated_cie_x);
    integral_Y = sum(spectral_array.*interpolated_cie_y);
    integral_Z = sum(spectral_array.*interpolated_cie_z);
    coordinate_x(i) = integral_X/(integral_X + integral_Y + integral_Z);
    coordinate_y(i) = integral_Y/(integral_X + integral_Y + integral_Z);
    coordinate_rgb(i,:) = xyz2rgb([coordinate_x(i) coordinate_y(i) (1 - coordinate_x(i) - coordinate_y(i))],'ColorSpace','adobe-rgb-1998');
    
    % Plot data. There are two options here - the first just plots the
    % curves normally, the second plots them in 3D with the third dimension
    % being the amount of strain for each spectrum. User should comment out the unused
    % one
    %plot(wavelength_array,100*spectral_array,'color',coordinate_rgb(i,:));
    plot3(wavelength_array,100*spectral_array,strain_array,'color',coordinate_rgb(i,:));   
end

% Set up axis ranges
if max(strain_array) == 0
    axis([wavelength_lower wavelength_upper 0 100 0 1]);
else
    axis([wavelength_lower wavelength_upper 0 100 0 max(strain_array)]);
end

% Set background color
set(gcf,'color','white');

% Set up CIE plot
figure(2);
hold on;
plotChromaticity;
set(gcf,'color','white');

% Create scatter plot, copy it to CIE plot, and close scatter plot
% This is a workaround, no other easy way to get scatter data onto CIE plot
figure(3);
hold on
plot(coordinate_x,coordinate_y,'w','LineWidth',1);
scatter(coordinate_x,coordinate_y,[],coordinate_rgb,'filled','LineWidth',1.5,'MarkerEdgeColor',[1 1 1]);
copyobj(findobj(3,'type','line'),findobj(2,'type','axes'));
copyobj(findobj(3,'type','scatter'),findobj(2,'type','axes'));
close(3);
