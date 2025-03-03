%% loop through the folder and organize the data.
% want it to be:
% run1: IMG2, IMG3, IMG4
% run2: IMG2, IMG3, IMG4
% ...
% and each image is 1040 by 1392
% 
% so it could be a 3D array for IMG2 (1040,1392,runs), 3D array for IMG3
% (1040, 1392, runs), 3D array for IMG4 (1040,1392,runs)
% and to make it all from the same run, call the same index of each (really
% will end up looping through the run #s and do the image processing in the
% loop

% Glass Cell: 1.56 um/pixel
% Main Chamber: 4.26 um/pixel

% IMG2: atoms and light
% IMG3: background- no atoms, no light
% IMG4: just the light

%% from Chat GPT
% Directory containing the files
folderPath = '/Users/Lauren/Documents/MATLAB/Side Image Analysis/20241211 data/images'; % Update this path

% List all CSV files in the folder
fileList = dir(fullfile(folderPath, '*.csv'));

% List all CSV files in the folder
%fileList = files;
disp(folderPath); % Display the folder path
fileList = dir(fullfile(folderPath, '*.csv'));
disp(fileList);   % Display the list of files


% Initialize containers for each IMG type
IMG2_data = [];
IMG3_data = [];
IMG4_data = [];

% Initialize a list to track unique run labels
runLabels = [];

% Process each file
for i = 1:length(fileList)
    % Get the filename
    fileName = fileList(i).name;

    % Parse the filename to extract IMG type and runlabel
    tokens = regexp(fileName, '^0-IMG(\d+)-([\d.]+)-.*\.csv$', 'tokens');
    if ~isempty(tokens)
        imgType = str2double(tokens{1}{1});     % IMG type (2, 3, or 4)
        runLabel = str2double(tokens{1}{2});   % Runlabel (e.g., 2.520000)

        % Add the runlabel to the list if it's new
        if ~ismember(runLabel, runLabels)
            runLabels = [runLabels, runLabel];
        end

        % Find the index of the runlabel in the sorted list
        [~, runIndex] = ismember(runLabel, sort(runLabels));

        % Read the data from the CSV file
        filePath = fullfile(folderPath, fileName);
        data = readmatrix(filePath); % Adjust if necessary

        % Store the data in the appropriate 3D array
        switch imgType
            case 2
                IMG2_data(:,:,runIndex) = data; % Assign to 3D array
            case 3
                IMG3_data(:,:,runIndex) = data;
            case 4
                IMG4_data(:,:,runIndex) = data;
        end
    else
        warning('Filename %s does not match the expected pattern.', fileName);
    end
end

%% Image processing
% start by just doing one set of images; later add a loop to loop through
% by run label.
% IMG2: atoms and light
% IMG3: background- no atoms, no light
% IMG4: just the light
od = zeros(1040,1392,i);
i = 1;
img2 = IMG2_data(:,:,i);
img3 = IMG3_data(:,:,i);
img4 = IMG4_data(:,:,i);
backsub4 = img4-img3;
backsub2 = img2-img3;
div42 = backsub4./backsub2;
od(:,:,i)=real(log(div42));
od1=od(:,:,i);

%% Calculating atom number
% now that we have optical depth, can estimate the atom number.
% following the LabVIEW script
sigmax = 0.346883; % scattering cross section of cesium in um^2
lambda = 0.852;
sig0 = 3*lambda^2/(2*3.141592);
I_sat = 1.1; % I_sat in mW per cm^2
linterm = (backsub4-backsub2)/I_sat;
% Glass Cell: 1.56 um/pixel
% Main Chamber: 4.26 um/pixel
pxconv = [4.26;1.56];

% from Chen Lung thesis: OD = nσ/(1 + M0/Msat)
% n = -(1/sigma)* ln(Pt/P0)
% where Pt is with atoms and light (backsub 2)
% and P0 is with just the light (backsub 4)
n1 = (1/sig0).*od1;
n2 = n1.*pxconv(1)^2;
% n2 is atoms per pixel.

%% need to set up a ROI and then fit a Gaussian to n2 in the roi
xroi = [250 1040]; % columns (1392), 2nd index (length 791)
yroi = [300 750]; % columns (1040), 1st index (length 451)
% when we plot, x goes across, y goes down
n2roi = n2(yroi(1):yroi(2),xroi(1):xroi(2));
% ok that worked. now integrate along each direction

n2xint = sum(n2roi,1); % to get x cut, sum along y
n2yint = sum(n2roi,2); % to get y cut, sum along x
x0 = 1:(xroi(2)-xroi(1))+1;
y0 = 1:(yroi(2)-yroi(1))+1;

% look at the integral linecuts
%%
figure;
t = tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
axis equal;
nexttile;
plot(x0, n2xint);
nexttile;
plot(y0,n2yint);

% now fit Gaussians of the form a0 + a1*x + a2*exp(-2*((x-a3)/a4)^2)

%% Dataset 1 % this is the x cross section
x1 = x0;
y1 = n2xint;

% Dataset 2 % this is the y cross section
x2 = y0;
y2 = n2yint;

% Define Gaussian fit type
gaussianModel = fittype('a0 + a1*x + a2*exp(-2*((x-a3)/a4)^2)', ...
                        'independent', 'x', ...
                        'coefficients', {'a0', 'a1', 'a2','a3','a4'});

% ('a*exp(-((x-b)^2)/(2*c^2))', ...
%                         'independent', 'x', ...
%                         'coefficients', {'a', 'b', 'c'});
% so a is like a2, b is like a3, c is like a4. need different guesses for a0 and a1


% Fit dataset 1 (the x axis data set)
fitResult1 = fit(x1(:), y1(:), gaussianModel, ...
                 'StartPoint', [min(y1),min(y1),max(y1), mean(x1), std(x1)]);

% Fit dataset 2 (the y axis data set)
fitResult2 = fit(x2(:), y2(:), gaussianModel, ...
                 'StartPoint', [min(y2),min(y2),max(y2), mean(x2), std(x2)]);

% Calculate FWHM
fwhmx = 2 * sqrt(2 * log(2)) * fitResult1.a4;
fwhmy = 2 * sqrt(2 * log(2)) * fitResult2.a4;

% Display results
disp('Dataset 1 Fit Parameters:');
disp(fitResult1);
disp('Dataset 2 Fit Parameters:');
disp(fitResult2);

% Plot the fits
figure;
subplot(2, 1, 1);
plot(fitResult1, x1, y1);
title('Gaussian Fit for Dataset 1');

subplot(2, 1, 2);
plot(fitResult2, x2, y2);
title('Gaussian Fit for Dataset 2');


%% Display IMG2, IMG3, IMG4, OD from Matlab, OD from LabVIEW
filePath = "/Users/Lauren/Documents/MATLAB/Side Image Analysis/20241211 data/images/a-od-2.520000-3.200000--0.900000.csv"
odlabview = readmatrix(filePath); % Adjust if necessary
%%
figure;
t = tiledlayout(1, 8, 'TileSpacing', 'Compact', 'Padding', 'Compact');
axis equal;
nexttile;
imagesc(img2);
title('IMG2')

nexttile;
imagesc(img3);
title('IMG3')

nexttile;
imagesc(img4);
title('IMG4')

nexttile;
imagesc(od(:,:,1));
title('OD Matlab')
colorbar;
clim([-0.5 1]);

nexttile;
imagesc(odlabview(:,:,1));
title('OD LabVIEW')
colorbar;
clim([-0.5 1]);


nexttile;
imagesc(n1);
title('n1')
colorbar;
clim([-0.5 2]);
axis equal
xlim(xroi);
ylim(yroi);

nexttile;
imagesc(n2);
title('n2')
colorbar;
axis equal;
clim([-0.5 2]);
xlim(xroi);
ylim(yroi);

nexttile;
imagesc(n2roi);
title('n2 roi')
colorbar;
axis equal;
clim([-0.5 2]);



