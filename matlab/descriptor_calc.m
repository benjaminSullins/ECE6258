function [ imageOrientation, imageEccentricity, imageWidth, imageLength, imageFingers, imageKnuckles,Fourier_mean,Fourier_max,Fourier_sigma,Fourier_min,Fourier_dc,Fourier_first] = descriptor_calc( imageStack )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

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
numFiles = size(imageStack,3);
imageStackMoments = zeros( (numMoments+1)^2, numFiles );


for i = 1:numFiles
    
    sz = size( imageStack(:,:,i) );
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
            Mpq = sum( reshape( bsxfun( @times, bsxfun( @times, imageStack(:,:,i), x.^p ), y.^q ), [], 1 ) );
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



%% Finger Segmentation
% Here, we find the contour outline of the hand and describe that as our
% boundary. We locate the wrist and measure the distances from the wrist to
% the boundary locations of the hand. Plotting this allows us to see where
% the finger locations are and help describe what the hand looks like.
% Analysis has shown that it might be better to use the number of
% fingertips and the number of knuckles. Fingertips have high standard
% deviations from the mean while the knuckles will be closer to the mean.
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
    
    % Attempting to find the fingertips through a search method across the
    % boundary for large curvature spikes.
    % Reference:
    % https://www.mathworks.com/matlabcentral/newsreader/view_thread/152405
    %
%     offset = 12;
%     threshold = 0.5;
%     
%     for l = 1 + offset:max(size(x)) - offset
%         xbounds = [1 2 3]; 
%         ybounds = [distance(l - 1) distance(l) distance(l + 1)];
%         
%         mx 	= mean(xbounds); 
%         my  = mean(ybounds);
%         
%         X = xbounds - mx; 
%         Y = ybounds - my; % Get differences from means
%         
%         dx2 = mean(X.^2); 
%         dy2 = mean(Y.^2); % Get variances
%         
%         t = [X,Y]\(X.^2-dx2+Y.^2-dy2)/2;    % Solve least mean squares problem
%         a0 = t(1); b0 = t(2);               % t is the 2 x 1 solution array [a0;b0]
%         r = sqrt(dx2+dy2+a0^2+b0^2);        % Calculate the radius
%         a = a0 + mx; b = b0 + my;           % Locate the circle's center
%         curv(l - offset) = 1/r;                         % Get the curvature
%         
%         % Threshold Curvature
%         if curv < threshold
%             peakX = [peakX x(l)];
%             peakY = [peakY y(l)];
%         end
%         
%     end

% Additional Fingertip classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% https://www.mathworks.com/matlabcentral/fileexchange/32696-2d-line-curvature-and-normals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    addpath( './linecurvature_version1b/' );
%    
%    sampling = 6;
%    
%    % x = (1:max(size(distance)))';
%    % y = distance;
%    
%    k = LineCurvature2D( [x(1:sampling:end) y(1:sampling:end)] );
%    k = upsample(k,sampling);
%    k = abs( k );
%    [ I, J, V ] = find( k > 0.09 );
%    figure;plot(k);
%    
%    N = LineNormals2D( [x y] );
%    
%    figure, imshow(binaryImage, []);
%    % figure, plot(x, y);
%    hold on;
%    plot(x(I), y(I), 'r^');
%    plot([x(1:sampling:end) x(1:sampling:end)+10*N(1:sampling:end,1)]',[y(1:sampling:end) y(1:sampling:end)+10*N(1:sampling:end,2)]');
%    hold off;

    % Now we find the maximum and minimum positions. These will be used to
    % determine if the finger is extended or not.
    [ peak peakIndexes ]        = findpeaks(distance, 'MinPeakHeight', 60);
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

% %     % Plotting assistance for viewing the finger distances.
%     figure; plot(distance);
%     hold on;
%     plot(peakIndexes, peak, 'r^');
%     plot(valleyIndexes, valley, 'go')
%     hold off;
%     
%     figure;
%     plot(curv);

% %     Plotting assistance for viewing the finger locations on the figure.
% %     figure, imshow(binaryImage, []);
% %     hold on;
% %     plot(x(peakIndexes), y(peakIndexes), 'r^');
% %     plot(x(valleyIndexes), y(valleyIndexes), 'go')
% %     hold off;
% % 
% %     This code can be used for construction the hull of the image. Think of it
% %     as the boundaries of the image only using the maximum difference
% %     locations. You can also plot the image if it helps make sense
% %     indexes = convhull(x, y);
% %     hold on; 
% %     plot(x(indexes), y(indexes), 'm-', 'LineWidth', 2);
% %     for k = 1 : length(peakIndexes)
% %         line([x(peakIndexes(k)), wristX], [y(peakIndexes(k)),wristY ], 'Color', 'r', 'LineWidth', 2);
% %     end
% %     hold off;
    
    % Store the Finger points as a descriptor
    imageFingers{1}{z} = max(size(peak));
    %imageFingers{2}{z} = peak;
    %imageFingers{3}{z} = peakIndexes;
    
    imageKnuckles{1}{z} = max(size(valley));
    %imageKnuckles{2}{z} = valley;
    %imageKnuckles{3}{z} = valleyIndexes;
    
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
tic
% Fourier Descriptors for the Images
Fourier_mean = zeros(1,numFiles);
Fourier_max = zeros(1,numFiles);
Fourier_sigma = zeros(1,numFiles);
Fourier_min = zeros(1,numFiles);
Fourier_dc =zeros(1,numFiles);
Fourier_first = zeros(1,numFiles);
for i = 1:numFiles
    image = imageStack(:,:,i);
    binaryImage = edge( image );
    [ avg, max_coeff ,sigma,min1,dc,firstharmonic ] = fourier_desc( binaryImage );
    Fourier_mean(i) = avg;
    Fourier_max(i) = max_coeff;
    Fourier_sigma(i) = sigma;
    Fourier_min(i) = min1;
    Fourier_dc(i) =dc;
    Fourier_first(i) = firstharmonic;
    


end
toc

