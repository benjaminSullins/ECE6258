function [ sign,corrPlot ] = classify( numVectors, trainingFeatureVector, C, varargin )
%This file computes the final sign output based on the training vector feature
%vectors and their corresponding group locations generated from the K-Means
%algorithm. A correlation is run amongst the estimated group regions to
%identify the most likely sign language character that the input belongs to.

display('Classifying Sign');
tic

numFiles = max(size(trainingFeatureVector(1,1,:)));

% Compute the feature vector for the input image.
% This is were we begin classifying the images amongst the training images
featureVector = zeros(1,numVectors);

for i = 1:numVectors
    
    % Convert if the underlying object is a cell
    if iscell( varargin{1}{i*2-1}(1) )
        X = cell2mat( varargin{1}{i*2-1}(1) );
        Y = cell2mat( varargin{1}{i*2}(1) );
    else
       X = varargin{1}{i*2-1}(1);
       Y = varargin{1}{i*2}(1); 
    end
    
    % Compute the distance from each centroid in the k-means clustering
    % algorithms. Here, X(1,:) represents the new sign.
    D = pdist2( [X Y], C{i}, 'euclidean');

    % The minimum distance represents the most likely group that the signal
    % belongs to.
    group = find( D == min(D) );

    % This group identifies what the key descriptor is tellings us what the
    % object is. This is added to the feature vector.
    featureVector(i) = group;

end

% Now that we have all of the feature vectors established for the input
% image, we can correlate this with the feature vectors of the trained
% images to see what we are most likely looking at.
% We determine what the sign is by correlating the group descriptors of the
% sign with the expected results of the training images. This reduces the
% overhead from a direct convolution while also trying to find the exact
% match.

currentCorrMax = 0;
sign = 9999;
corrPlot = zeros(1,numFiles);

for i = 1:numFiles
    
    correlation = xcorr( featureVector, trainingFeatureVector(:,:,i), 'coeff' );
    
    % Find the correlation value at the central location of the feature
    % vector array.
    correlation = correlation(numVectors);
    
    % For debugging purposes, we plot the correlation amongst the multitude
    % of signs.
    corrPlot(i) = correlation;
    
    % Check if the current sign is more likely to be the match
    if correlation > currentCorrMax
        currentCorrMax  = correlation;
        sign            = i;
    end
end

% % Plot the correlation
% figure;
% plot( corrPlot,'k*','MarkerSize',5);

% The input image has been found based on it's correlation with the
% training image feature vectors.
toc;

end

