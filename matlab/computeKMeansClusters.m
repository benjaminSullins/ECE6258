function [ C ] = computeKMeansClusters( numVectors, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trying to do K Means
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numVectors
    
    % Convert if the underlying object is a cell
    if iscell( varargin{1}{i*2-1} )
        varargin{1}{i*2-1} = cell2mat( varargin{1}{i*2-1} );
        varargin{1}{i*2} = cell2mat( varargin{1}{i*2} );
    end
    
    X = [ varargin{1}{i*2-1}' varargin{1}{i*2}'];

    % Find the predicted and actual centroid locations
    opts = statset('Display','final');
    [idx,centroids] = kmeans(X,5,'Distance','sqeuclidean',...
        'Replicates',5,'Options',opts);
    
    % Store the centroids
    C{i} = centroids;
    
    % Defines a fine grid on the plot
    
    step = 0;
    for p = 1 : -0.001 : 0
       if max(size( min(X(:,1)) : p : max(X(:,1)) )) > 500
           step = p;
           break;
       end
    end
    
    x1 = min(X(:,1)):step:max(X(:,1));
    x2 = min(X(:,2)):step:max(X(:,2));
    [x1G,x2G] = meshgrid(x1,x2);
    XGrid = [x1G(:),x2G(:)];

    idx2Region = kmeans(XGrid,5,'MaxIter',1,'Start',C{i});

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
end

end

