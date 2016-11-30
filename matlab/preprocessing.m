function [ imageStack,check ] = preprocessing( filePath, filename, key )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
%key = 0 Training phase, load multiple images
%key = 1 Testing phase, load single image
%% Import Training Images
disp('Importing Training Images');

% filePath = '../images/';
files = dir(filePath);
[numFiles null] = size(files);


% The images are not all the same size. We wil resize the images and then
% modify the algorithms to handle the discrepencies.
imagePaddSize = [200 200];

% Read out all the files. The first few inputs are not files, so we offset
% to number 3 instead.
offset = 3;
check = numFiles;

%Read single image
if key ==1
    if filename ~= 0
        numFiles = 1;
    else
        check =0;
    end
    
    offset =1;
    files(1).name = filename;
end

imageStack = zeros(imagePaddSize(1),imagePaddSize(2),numFiles-offset+1);
%if numFiles<offset
%    return;
    
for i = offset:numFiles

    % Read image
    image = imread(strcat(filePath, files(i).name), 'jpg');
    
    % Convert to grayscale
    image = rgb2gray(image);
    
    image = imresize(image, imagePaddSize);
    
    image = intmax(class(image))  - image; % We want to flip the scale of our images because they are black hot.
    
    image = im2double(image);
    
    % Add to image stack
    imageStack(:,:,i - (offset - 1)) = image;
    
    
end
end

