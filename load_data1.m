
clear; clc; close all

data_path = "data\2024\12 December\20241211\0333-";
files = dir(data_path + "/images/*.csv");

%%
data = readmatrix(['data\2024\12 December\20241211\0333-', '/images/', '0-IMG4-2.420000-3.200000--0.900000.csv']);
imagesc2(data)

clear; clc; close all

data_path = "/Users/Lauren/Documents/MATLAB/Side Image Analysis/20241211 data";
files = dir(data_path + "/images/*.csv");

%%
data = readmatrix(['20241211 data/images/', '0-IMG2-2.520000-3.200000--0.900000.csv']);
imagesc(data)

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

