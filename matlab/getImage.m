function [ image fileName ] = getImage( filePath, fileNumber )
%Imports and image and returns it to the calling function.
%   Detailed explanation goes here

%% Getting Image
disp('Importing Image');

% filePath = '../images/';
files = dir(filePath);
[numFiles null] = size(files);


% The images are not all the same size. We wil resize the images and then
% modify the algorithms to handle the discrepencies.
imagePaddSize = [200 200];

% Read out all the files. The first few inputs are not files, so we offset
% to number 2 instead.
offset = 2;

% Read image
image = imread(strcat(filePath, files(fileNumber + offset).name), 'jpg');

fileName = files(fileNumber + offset).name;

end

