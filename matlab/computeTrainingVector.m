function [ trainingFeatureVector ] = computeTrainingVector( numVectors, C, varargin )
%This file takes the feature vector extractions and generates an array for
%future execution amongst the results. The file also computes the distance
%between the K-Means centroids to selectively choose which centroid the
%training vector belongs to.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the training feature vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numFiles = max(size(varargin{1}{1}));

% Generate our training vector list
trainingVector = zeros( numFiles, numVectors*2);

for i = 1:numFiles
    trainingVector(i, 1) = varargin{1}{1}(i); % imageOrientation(i);
    trainingVector(i, 2) = varargin{1}{2}(i); % imageEccentricity(i);
    trainingVector(i, 3) = varargin{1}{3}(i); % imageLength(i);
    trainingVector(i, 4) = varargin{1}{4}(i); % imageWidth(i);
    trainingVector(i, 5) = varargin{1}{5}(i);
    trainingVector(i, 6) = varargin{1}{6}(i);
    trainingVector(i, 7) = varargin{1}{7}(i);
    trainingVector(i, 8) = varargin{1}{8}(i);
    trainingVector(i, 9) = varargin{1}{9}(i);
    trainingVector(i, 10) = varargin{1}{10}(i);
    trainingVector(i, 11) = varargin{1}{11}(i);
    trainingVector(i, 12) = varargin{1}{12}(i);
%     trainingVector(i, 13) = varargin{1}{13}(i);
%     trainingVector(i, 14) = varargin{1}{14}(i);
%     trainingVector(i, 15) = varargin{1}{15}(i);
%     trainingVector(i, 16) = varargin{1}{16}(i);
end

% Compute the feature vectors for all of the training images
trainingFeatureVector = zeros(1,numVectors,numFiles);

for i = 1:numFiles
    % For each vector
   for j = 1:numVectors
      
       % Compute the distance from each centroid in the k-means clustering
       % algorithms. Here, X(1,:) represents the new sign.
       D = pdist2( [trainingVector(i,j*2-1) trainingVector(i,j*2)], C{j}, 'euclidean');

       % The minimum distance represents the most likely group that the signal
       % belongs to.
       group = find( D == min(D) );
       
       trainingFeatureVector(1,j,i) = group;
       
   end
end

end

