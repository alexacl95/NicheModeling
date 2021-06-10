function ReadInfo = ReadLayers(layerfolder,parallel)
    tic
    if nargin < 2
        parallel = false;
    end

    if parallel
        parforArg = inf;
    else
        parforArg = 0;
    end
    
    %Layer lecture
    addpath(layerfolder)
    
    LayerDir = dir(layerfolder);
    layers = {LayerDir.name};
    CompareFileLayers = strncmp('bio', layers, 3);
    layers = layers(CompareFileLayers);
    [Z,R] = arcgridread(strcat(layerfolder,layers{1}));
    N = length(layers);
    Z = zeros(size(Z, 1), size(Z, 2),N);

    disp('----Reading layers----')
    
    parfor (i = 1 : N, parforArg)
        %progressbar(i/N)
        interZ = readgeoraster(layers{i},'CoordinateSystemType','geographic');
        %interZ = arcgridread(layers{i});
        interZ(interZ==-9999)=NaN;
        Z(:, :, i)=interZ;
    end
    %Create a new data set 
    [Dimension1, Dimension2, Dimension3] = size(Z);
    Map = Z(:, :, 1);
    ClimVariables = zeros(Dimension3, Dimension1 * Dimension2);
    
    parfor (i = 1 : Dimension3, parforArg)
        aux = Z(:, :, i);
        ClimVariables(i, :) = aux(:);
    end
    
    rmpath(layerfolder)
    
    nanDetector = sum(ClimVariables);
    nanPositions = isnan(nanDetector);
    ClimVar = ClimVariables(:,~nanPositions);
    Dimension = size(ClimVar,2);
    Map(:) = mean(ClimVariables);
    
    ReadInfo.Indicator = nanPositions;
    ReadInfo.Dimensions = [Dimension, Dimension3];
    %ReadInfo.Distance = zeros(1, Dimension1 * Dimension2);
    ReadInfo.NormalizedClimVar = normalize(ClimVar, 2, 'range');
    ReadInfo.Map = Map;
    ReadInfo.Z = Z;
    ReadInfo.R = R;
    toc
    
end