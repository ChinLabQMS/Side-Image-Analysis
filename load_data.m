clear; clc; close all

data_path = "data\2024\12 December\20241211\0333-";
files = dir(data_path + "/images/*.csv");

%%
data = readmatrix(['data\2024\12 December\20241211\0333-', '/images/', '0-IMG4-2.420000-3.200000--0.900000.csv']);
imagesc2(data)
