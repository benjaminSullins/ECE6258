% ECE 6258: Digital image processing 
% Benjamin Sullins
% GTID: 903232988
% Distance Learning Student
% School of Electrical and Computer Engineering 
% Georgia Instiute of Technology 
%
% Final Project
% Sign Language Translation To Text
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date Modified : 11/22/2016 - Ben Sullins
% - Initial design implementation
% Date Modified : 11/24/2016 - Ben Sullins
% - Added Fingertip Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References
% http://cs229.stanford.edu/proj2011/KurdyumovHoNg-SignLanguageClassificationUsingWebcamImages.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
fclose('all');

%% Import Training Images
disp('Importing Training Images');
tic

filePath = '../images/';
files = dir(filePath);
[numFiles null] = size(files);

% The images are not all the same size. We wil resize the images and then
% modify the algorithms to handle the discrepencies.
imagePaddSize = [200 200];

% Read out all the files. The first few inputs are not files, so we offset
% to number 3 instead.
offset = 3;
for i = offset:numFiles

    % Read image
    image = imread(strcat(filePath, files(i).name), 'jpg');
    
    % Convert to grayscale
    image = rgb2gray(image);
    image = imresize(image, imagePaddSize);
    
    % Add to image stack
    imageStack(:,:,i - (offset - 1)) = image;
    
end

[imSizeY imSizeX numFiles] = size(imageStack);

clear filePath files null null2 imagePaddSize image i offset
toc

% We want to flip the scale of our images because they are black hot.
disp('Converting Images to White Hot');
tic
for i = 1:numFiles
   imageStack(:,:,i) = intmax(class(imageStack)) - imageStack(:,:,i);
end

clear i
toc

%% Intertial Moment Calculation
% The first set of key descriptors will be the moments (inertia)
% calculations of the image
%
% References:
% http://stackoverflow.com/questions/27684478/how-do-i-calculate-the-central-moment-in-an-image
% https://en.wikipedia.org/wiki/Image_moment
% https://en.wikipedia.org/wiki/Eccentricity_(mathematics)
% http://breckon.eu/toby/teaching/dip/opencv/SimpleImageAnalysisbyMoments.pdf
    
disp('Computing Moment Calculations');
tic

% Compute up to the third moments
numMoments = 3;
imageStackMoments = zeros( (numMoments+1)^2, numFiles );

% Convert the images to doubles for computation
for i = 1:numFiles
   imageStackDoubles(:,:,i) = im2double( imageStack(:,:,i) );
end

for i = 1:numFiles
    
    sz = size( imageStackDoubles(:,:,i) );
    x = ( 1:sz(2) );
    y = ( 1:sz(1) ).';
    x = x - mean(x);
    y = y - mean(y);

    % Compute the p-q moments
    for p = 0:numMoments
        for q = 0:numMoments 
            
            % Skipping Conditions
            if (p == 1 || p == 2) && q == 3
                % Skip (1,3), (2,3)
                continue;
            elseif p == 2 && q == 2
                % Skip the (2,2) moment
                continue;
            elseif p == 3
                % Break all after (3,0)
                if q ~= 0
                    break;
                end
            end
            
            % Compute the moments
            Mpq = sum( reshape( bsxfun( @times, bsxfun( @times, imageStackDoubles(:,:,i), x.^p ), y.^q ), [], 1 ) );
            imageStackMoments( p*(numMoments+1) + q + 1,i) = Mpq;
        end
    end
    
end

% Compute the Moments
% The moments can be used to find defining descriptors of the image.
% Use the following source for a further description :
% http://breckon.eu/toby/teaching/dip/opencv/SimpleImageAnalysisbyMoments.pdf

imageOrientation    = zeros(1, numFiles);
imageEccentricity   = zeros(1, numFiles);

for i = 1:numFiles
    % Convert the moments to the Central Moments. This removes the
    % translation dependency.
    % See page 3
    m10Prime = 0;
    m01Prime = 0;
	m20Prime = imageStackMoments(2*(numMoments+1)+0+1, i) / imageStackMoments(0*(numMoments+1)+0+1, i) ...
               - (imageStackMoments(1*(numMoments+1)+0+1, i) / imageStackMoments(0*(numMoments+1)+0+1, i))^2; 
    m02Prime = imageStackMoments(0*(numMoments+1)+2+1, i) / imageStackMoments(0*(numMoments+1)+0+1, i) ...
               - (imageStackMoments(0*(numMoments+1)+1+1, i) / imageStackMoments(0*(numMoments+1)+0+1, i))^2;  
    m11Prime = imageStackMoments(1*(numMoments+1)+1+1, i) / imageStackMoments(0*(numMoments+1)+0+1, i) ...
               - (imageStackMoments(1*(numMoments+1)+0+1, i) / imageStackMoments(0*(numMoments+1)+0+1, i)) ...
               * (imageStackMoments(0*(numMoments+1)+1+1, i) / imageStackMoments(0*(numMoments+1)+0+1, i));
    
    % Compute the Orientation (Direction of the image)
    % The orientation of the object is defined as the tilt angle between the x-axes and
    % the axis, around which the object can be rotated with minimal inertia
    % Check zero condition
    if m20Prime ~= m02Prime
        imageOrientation(i) = (0.5) * atan( (2 * m11Prime) / (m20Prime - m02Prime) );
    end
    
    % Compute the Eccentricity (How elongated the image is)
    % The eccentricity ? can have values from 0 to 1. It is 0 with a perfectly round
    % object and 1 by a line shaped object.
    
    % The eccentricity ? is a better measure than the roundness ? of the object,
    % because it has a clearly defined range of values and therefore it can be compared
    % much better.
    imageEccentricity(i) = ( ((m20Prime - m02Prime)^2) - 4*(m11Prime^2) ) / (m20Prime + m02Prime)^2;
    
end

clear sz x y p q Mpq i numMoments m20Prime m02Prime m11Prime lambda1 lambda2 ...
        m10Prime m01Prime imageStackDoubles
toc

%% Finger Segmentation
% Here, we find the contour outline of the hand and describe that as our
% boundary. We locate the wrist and measure the distances from the wrist to
% the boundary locations of the hand. Plotting this allows us to see where
% the finger locations are and help describe what the hand looks like.
% Analysis has shown that it might be better to use the fingertip peaks and
% not the valleys due to the errors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% http://www.mathworks.com/matlabcentral/answers/215225-can-we-plot-the-boundary-of-an-image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing Finger Segmentations');
tic

for z = 1:numFiles

    % Find the contour of the image.
    % This may need to be modified to a true countour instead of thresholding
    % with the im2bw() function.
    binaryImage = im2bw( imageStack(:,:,z) );
    % binaryImage = contour(binaryImage);

    % figure, imshow(binaryImage, []);

    % Find the boundary locations of the countoured image
    boundaries = bwboundaries(binaryImage, 'noholes');
    x = boundaries{1}(:, 2);
    y = boundaries{1}(:, 1);

    % We now have the boundary locations for the contour. The wrist is usually
    % located at the bottom of the image, so we find the bottom most y location
    % in the contour. The x location at the max y locations is found and
    % averaged to help place the central point in the middle of the hand's
    % wrist.
    % This could be improved by finding the smallest difference in x posiitons
    % as this would locate the wrist on most people.
    wristY = max(y);
    wristX = ( max(x(find( y == max(y) ))) + min(x(find( y == max(y) ))) )/2;

    % Compute the distances from the central wrist location to the boundary
    % locations along the contour.
    distance = sqrt( (x - wristX).^2 + (y - wristY).^2 );

    % Now we find the maximum and minimum positions. These will be used to
    % determine if the finger is extended or not.
    [ peak peakIndexes ]        = findpeaks(distance);
    [ valley valleyIndexes ]    = findpeaks(max(distance) - distance);
    valley = max(valley) - valley;

    % The results will contain some erroneous data as it's sensitive to changes
    % in the distance. Rules will be applied and filtering techniques to help
    % remove those erroneous measures.
    % Rule 1: For most cases, the finger tips are the farthest from the wrist.
    % This allows us to set a threshold to limit any peaks/valleys which are
    % close to the wrist location.
    threshold       = 60;
    peakIndexes     = peakIndexes(find( peak > threshold ));
    peak            = peak(find( peak > threshold ));

    valleyIndexes   = valleyIndexes(find( valley > threshold ));
    valley          = valley( valley > threshold );

    % Rule 2: Fingertips are usually spaced out with some distance between
    % them. When there are multiple peaks next to each other, the maximum one
    % is selected as the finger location.
    threshold = 10;

    % Cycle through all the peaks
    temp = peak;
    for i = 1:size(peak)-1

        % If the two peak values are close to each other
        if abs(peak(i+1) - peak(i)) < threshold
            % Find the largest peak
            if peak(i+1) > peak(i)
                % Remove the smaller
                temp(i) = 0;
            else
                % Remove the smaller
                temp(i+1) = 0;
            end
        end
    end
    peak = temp;

    % Remove the cleared peaks
    peakIndexes = peakIndexes(find( peak > 0 ));
    peak        = peak(find( peak > 0 ));

    % Cycle through all the valleys
    temp = valley;
    for i = 1:size(valley)-1
        % If the two valley values are close to each other
        if abs(valley(i+1) - valley(i)) < threshold
            % Find the smallest valley
            if valley(i+1) < valley(i)
                % Remove the larger
                temp(i) = 0;
            else
                % Remove the larger
                temp(i+1) = 0;
            end
        end
    end
    valley = temp;

    % Remove the cleared peaks
    valleyIndexes   = valleyIndexes(find( valley > 0 ));
    valley          = valley(find( valley > 0 ));

    % Plotting assistance for viewing the finger distances.
%     figure; plot(distance);
%     hold on;
%     plot(peakIndexes, peak, 'r^');
%     plot(valleyIndexes, valley, 'go')
%     hold off;

%     % Plotting assistance for viewing the finger locations on the figure.
%     figure, imshow(binaryImage, []);
%     hold on;
%     plot(x(peakIndexes), y(peakIndexes), 'r^');
%     plot(x(valleyIndexes), y(valleyIndexes), 'go')
%     hold off;
% 
%     % This code can be used for construction the hull of the image. Think of it
%     % as the boundaries of the image only using the maximum difference
%     % locations. You can also plot the image if it helps make sense
%     indexes = convhull(x, y);
%     hold on; 
%     plot(x(indexes), y(indexes), 'm-', 'LineWidth', 2);
%     for k = 1 : length(peakIndexes)
%         line([x(peakIndexes(k)), wristX], [y(peakIndexes(k)),wristY ], 'Color', 'r', 'LineWidth', 2);
%     end
%     hold off;
    
    % Store the Finger points as a descriptor
    imageFingers{1}{z} = size(peak);    % How many fingers are there?
    imageFingers{2}{z} = peak;          % Those finger values
    imageFingers{3}{z} = peakIndexes;   % Those finger positions
    
end

clear boundaries distance i indexes k peak peakIndexes temp threshold ...
        valley valleyIndexes wristX wristY x y z binaryImage

toc

%% Length and Width
% Compute the length and width of the hand for another key descriptor
disp('Computing Length and Width');
tic

imageWidth  = zeros(1, numFiles);
imageLength = zeros(1, numFiles);

for z = 1:numFiles
    
    % Find the contour of the image.
    % This may need to be modified to a true countour instead of thresholding
    % with the im2bw() function.
    binaryImage = im2bw( imageStack(:,:,z) );
    % binaryImage = contour(binaryImage);

    % Find the boundary locations of the countoured image
    boundaries = bwboundaries(binaryImage, 'noholes');
    x = boundaries{1}(:, 2);
    y = boundaries{1}(:, 1);
    
    % Compute the hull of the contour
    indexes = convhull(x, y);
    
    imageWidth(z) = (max(x(indexes)) - min(x(indexes)) );
    imageLength(z) = (max(y(indexes)) - min(y(indexes)) );
    
    midX = (max(x(indexes)) + min(x(indexes)) )/2;
    midY = (max(y(indexes)) + min(y(indexes)) )/2;
    
    % This code can be used for construction the hull of the image. Think of it
    % as the boundaries of the image only using the maximum difference
    % locations. You can also plot the image if it helps make sense
%     figure, imshow(binaryImage, []);
%     hold on; 
%     plot(x(indexes), y(indexes), 'm-', 'LineWidth', 2);
%     for k = 1 : length(imageWidth)
%         line([max(x(indexes)), min(x(indexes))], [midY, midY ], 'Color', 'r', 'LineWidth', 2);
%         line([midX, midX ], [max(y(indexes)), min(y(indexes))], 'Color', 'r', 'LineWidth', 2);
%     end
%     hold off;
    
end

clear x y z indexes boundaries binaryImage midX midY

toc

%% K Means Clustering Algorithm
% With all of the key descriptors calculated, the keys can be used to build
% a function describing the discrepencies between the signs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% https://www.mathworks.com/help/stats/kmeans.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
disp('K Means Clustering');
tic

% Plot the Orientation vs the Eccentricity
% figure;
% plot( imageOrientation, imageEccentricity,'k*','MarkerSize',5);

% Plot the Length vs. Width
% figure;
% plot( imageWidth, imageLength,'k*','MarkerSize',5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trying to do K Means on Length and Width
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = [ imageWidth' imageLength'];

% Find the predicted and actual centroid locations
opts = statset('Display','final');
[idx,C] = kmeans(X,5,'Distance','sqeuclidean',...
    'Replicates',5,'Options',opts);

% Defines a fine grid on the plot
x1 = min(X(:,1)):0.1:max(X(:,1));
x2 = min(X(:,2)):0.1:max(X(:,2));
[x1G,x2G] = meshgrid(x1,x2);
XGrid = [x1G(:),x2G(:)];

idx2Region = kmeans(XGrid,5,'MaxIter',1,'Start',C);

% Assigns each node in the grid to the closest centroid    
figure;
gscatter(XGrid(:,1),XGrid(:,2),idx2Region,...
    [0,0.75,0.75;0.75,0,0.75;0.75,0.75,0],'..');
hold on;
plot(X(:,1),X(:,2),'k*','MarkerSize',5);
title 'Width vs. Length';
xlabel 'Width (pixels)';
ylabel 'Length (pixels)';
legend('Region 1','Region 2','Region 3','Region 4','Region 5','Data','Location','SouthEast');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trying to do K Means on Orientation and Eccentricity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = [ imageOrientation' imageEccentricity' ];

% Find the predicted and actual centroid locations
opts = statset('Display','final');
[idx,C] = kmeans(X,5,'Distance','sqeuclidean',...
    'Replicates',5,'Options',opts);

% Defines a fine grid on the plot
x1 = min(X(:,1)):0.001:max(X(:,1));
x2 = min(X(:,2)):0.001:max(X(:,2));
[x1G,x2G] = meshgrid(x1,x2);
XGrid = [x1G(:),x2G(:)];

idx2Region = kmeans(XGrid,5,'MaxIter',1,'Start',C);

% Assigns each node in the grid to the closest centroid    
figure;
gscatter(XGrid(:,1),XGrid(:,2),idx2Region,...
    [0,0.75,0.75;0.75,0,0.75;0.75,0.75,0],'..');
hold on;
plot(X(:,1),X(:,2),'k*','MarkerSize',5);
title 'Orientation vs. Eccentricity';
xlabel 'Orientation';
ylabel 'Eccentricity';
legend('Region 1','Region 2','Region 3','Region 4','Region 5','Data','Location','SouthEast');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If we are given an input, how do we classify the image based on our key
% Descriptors?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate our training vector list
numVectors = 2;
trainingVector = zeros(numFiles, numVectors*2);

for i = 1:numFiles
    trainingVector(i, 1) = imageOrientation(i);
    trainingVector(i, 2) = imageEccentricity(i);
    trainingVector(i, 3) = imageLength(i);
    trainingVector(i, 4) = imageWidth(i);
end

% Compute the feature vectors for all of the training images
trainingFeatureVector = zeros(1,numVectors,numFiles);

for i = 1:numFiles
    % For each vector
   for j = 1:numVectors
      
       % Compute the distance from each centroid in the k-means clustering
       % algorithms. Here, X(1,:) represents the new sign.
       D = pdist2( [trainingVector(i,j*2-1) trainingVector(i,j*2)], C, 'euclidean');

       % The minimum distance represents the most likely group that the signal
       % belongs to.
       group = find( D == min(D) );
       
       trainingFeatureVector(1,j,i) = group;
       
   end
end

% Compute the feature vector for the input image.
% This is were we begin classifying the images amongst the training images
featureVector = zeros(1,numVectors);

for i = 1:numVectors
    
    % Compute the distance from each centroid in the k-means clustering
    % algorithms. Here, X(1,:) represents the new sign.
    D = pdist2(X(1,:), C, 'euclidean');

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

% Plot the correlation
figure;
plot( corrPlot,'k*','MarkerSize',5);

% The input image has been found based on it's correlation with the
% training image feature vectors.
toc;