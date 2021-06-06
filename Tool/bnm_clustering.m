function T = bnm_clustering(T, opt)

    %Species classification using spectral clustering
    
    %Setting the clustering options
    if nargin < 2        
        DiffMag = 1;
        Distance = "mahalanobis";
        InitialCluster = 10;
        ClassVar = 4:22;
        
    else        
        DiffMag = opt.DiffMag;
        Distance = opt.Distance;
        InitialCluster = opt.InitialCluster;
        ClassVar = opt.ClassVar;
        
    end

    X = table2array(T(:, ClassVar));
    [~, ~, D_temp] = spectralcluster(X, InitialCluster, 'Distance', Distance);
    
    % find order of magnitude
    n = floor(log(abs(D_temp))./log(10));
    INF = isinf(n);
    n(INF) = 0;
    
    % find optmial number of clusters
    k = find(n <= min(n) + DiffMag);
    k = length(k) + sum(INF);
    
    if ~isempty(k)        
        %Clasiffy each point into a cluster
        [ClusterIndex, ~, ~] = spectralcluster(X, k,'Distance', Distance);        
    else         
       k = 1;
       ClusterIndex = ones(1, length(T.LAT));
       disp("Warning: check the spectralcluster options")
    end
    
    %Save the classification index for each point
    %T = addvars(T, ClusterIndex,'Before','LAT');
    T = addprop(T, {'ClusterIndex', 'NumClusters'}, {'table', 'table'});
    T.Properties.CustomProperties.ClusterIndex = ClusterIndex;
    T.Properties.CustomProperties.NumClusters = k;
    
end