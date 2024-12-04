% Script to read in the csv file of the side CCD images
% Lauren Weiss, August 2024
% import data
% data = load('C:\Users\qmspc\Desktop\NewLabData\2024\08 August\test2.csv');
% 
% data2 = load('C:\Users\qmspc\Desktop\NewLabData\2024\08 August\test3.csv');
% 
% data5 = load('C:\Users\qmspc\Desktop\NewLabData\2024\08 August\test5.csv');

data9 = load('C:\Users\qmspc\Desktop\NewLabData\2024\08 August\test9.csv');

%% 
datado = data9;
clims = [-0.5 1]
figure
imagesc(datado,clims)
colorbar()

%% 
% unpack the first data set into its constituent images
numImg = numel(datado)/1392/1040;
img = zeros(1040, 1392, numImg);
for i = 1:numImg
    img(:,:,i) = datado((i*1040-1039):1040*i,:);
end
clims = [-0.5 1]


figure
imagesc(img(:,:,1),clims)
colorbar()

% so it seems like it's ok to pack the images together in this way.See if I
% can create a blank file in LabVIEW so I don't have to write them to a
% specific file...

%% 
% unpack the desired next data set into its constituent images
datado = data9;
numImg = numel(datado)/1392/1040;
img = zeros(1040, 1392, numImg);
for i = 1:numImg
    img(:,:,i) = datado((i*1040-1039):1040*i,:);
end
clims = [-0.5 1]

%%
figure
imagesc(img(:,:,2),clims)
colorbar()

%% 
clims = [-0.5 1]
figure
imagesc(datado,clims)
colorbar()

%% Reading binary data file... or at least trying to. it's not working
binarydata ='C:\Users\qmspc\Desktop\NewLabData\2024\08 August\test_export2.csv';

fileID = fopen(binarydata);
A = fread(fileID, 'uint32');
fclose(fileID);


dim1 = round(length(A)/1392);
a = zeros(round(length(A)/1392),1392);
for i = 1:round(length(A)/1392)
    a(i,:)=A(i*1392-1391:i*1392);
end

datado = a;
numImg = numel(datado)/1392/1040;
img = zeros(1040, 1392, numImg);
for i = 1:numImg
    img(:,:,i) = datado((i*1040-1039):1040*i,:);
end
clims = [-0.5 1]

%%
figure
imagesc(datado)
colorbar()

%%
databintry = labviewload('C:\Users\qmspc\Desktop\NewLabData\2024\08 August\test_export4.csv',2)
clims = [-0.5 1]
figure
imagesc(cell2mat(databintry),clims)
%% 
% LABVIEWLOAD  Load Labview binary data into Matlab
%
%   DATA = LABVIEWLOAD(FNAME,DIM); % Loads the Labview data in binary
%      format from the file specified by FNAME, given the NUMBER of 
%      dimensions (not the actual dimensions of the data in the binary 
%      file) of dimensions specified by DIM.
%
%      LABVIEWLOAD will repeatedly grab data until the end of the file.
%      Labview arrays of the same number of dimensions can be repeatedly 
%      appended to the same binary file.  Labview arrays of any dimensions
%      can be read.
%
%   DATA = LABVIEWLOAD(FNAME,DIM,PREC); % Loads the data with the specified
%      precision, PREC.
%
%   Note:  This script assumes the default parameters were used in the
%      Labview VI "Write to Binary File".  Labview uses the Big Endian
%      binary number format.
%
%   Examples:
%       D = labviewload('Data.bin',2);  % Loads in Data.bin assuming it
%          contains double precision data and two dimensions.
%
%       D = labviewload('OthereData.bin',3,'int8');  % Loads in
%          OtherData.bin assuming it contains 8 bit integer values or
%          boolean values.
%
% Jeremiah Smith
% 4/8/10
% Last Edit: 5/6/10


function data = labviewload(fname,dim,varargin)
siz = [2^32 2^16 2^8 1]';  % Array dimension conversion table

% Precision Input
if nargin == 2
    prec = 'double';
elseif nargin == 3
    prec = varargin{1};
else
    error('Too many inputs.')
end

%% Initialize Values

fid = fopen(fname,'r','ieee-be');  % Open for reading and set to big-endian binary format
fsize = dir(fname);  % File information
fsize = fsize.bytes;  % Files size in bytes

%% Preallocation
rows = [];
columns = [];
I = 0;
while fsize ~= ftell(fid)
    dims = [];
    for i=1:1:dim
        temp = fread(fid,4);
        temp = sum(siz.*temp);
        dims = [dims,temp];
    end
    I = I + 1;
%     fseek(fid,prod(dims)*8,'cof');  % Skip the actual data
    temp = fread(fid,prod(dims),prec,0,'ieee-be');  % Skip the actual data (much faster for some reason)
end
fseek(fid,0,'bof'); % Reset the cursor
data = repmat({NaN*ones(dims)},I,1);  % Preallocate space, assumes each section is the same

%% Load and parse data
for j=1:1:I
    dims = [];  % Actual array dimensions
    for i=1:1:dim
        temp = fread(fid,4);
        temp = sum(siz.*temp);
        dims = [dims,temp];
    end
    clear temp i
    
    temp = fread(fid,prod(dims),prec,0,'ieee-be');  % 11 is the values per row,
        % double is the data type, 0 is the bytes to skip, and
        % ieee-be specified big endian binary data
    
%% Reshape the data into the correct array configuration
    if dim == 1
        temp = reshape(temp,1,dims);
    else
        evalfunc = 'temp = reshape(temp';
        for i=1:1:dim
            evalfunc = [evalfunc ',' num2str(dims(dim-i+1))];
        end
        if dim ~= 2
            eval([evalfunc ');'])
        else
            eval([evalfunc ')'';'])
        end
    end
    
    data{j} = temp;  % Save the data
end

fclose(fid);  % Close the file
end