clear; clc; close all

data_path = "data/2024/12 December/20241203/0837-MOT loading";
files = dir(data_path + "/images/*.csv");

Data = struct('IMG2', [], 'IMG3', [], 'IMG4', []);

for i = 1: length(files)
    f = files(i);
    if contains(f.name, 'IMG2')
        category = "IMG2";
    elseif contains(f.name, 'IMG3')
        category = "IMG3";
    elseif contains(f.name, 'IMG4')
        category = "IMG4";
    elseif contains(f.name, 'OD')
        category = "OD";
    else
        warning('Unable to categorize image!')
    end
    if category == "OD"
    else
        params = sscanf(f.name, "%d-" + category + "-%f-%f-%f.csv");
        shot = struct('data', readmatrix([f.folder, '/', f.name]), 'scan1', params(2), ...
            'scan2', params(1), 'scan3', params(3), 'scan4', params(4));
        Data.(category) = [Data.(category), shot];
    end
end

%%
imagesc2(Data.IMG2(1).data)
