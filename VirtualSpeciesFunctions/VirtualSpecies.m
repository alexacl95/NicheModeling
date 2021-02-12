function [T, Map] = VirtualSpecies(ReadInfo, InfoInitialPoint, samples, Occupation, factor, show, spName)
        
    %%%-- Reading clim variables--%%

    Indicator = ReadInfo.Indicator;
    Dimension = ReadInfo.Dimensions(1);
    R = ReadInfo.R;
    [Dim1, Dim2] = size(ReadInfo.Map);
    Map = nan(Dim1, Dim2);
        
    NormDistance = InfoInitialPoint.NormDistance;
    idx = InfoInitialPoint.idx;
    SortNormDistance = InfoInitialPoint.SortNormDistance;
    
    limit = round(Dimension * Occupation); 
    
    SortNormDistance(limit : end) = 0;
    IndexSortNorm = 1 : limit - 1;
    SortNormDistance(IndexSortNorm) = normalize(SortNormDistance(IndexSortNorm), 2, 'range');
    
    NormDistance(idx) = SortNormDistance;
    Map(~Indicator) = NormDistance;
    
%     if show == 1
%         clf
%         geoshow(Map, R, 'DisplayType', 'surface');
%         contourcmap('jet', 0 : 0.05 : 1, 'colorbar', 'on', 'location', 'vertical')
%     end
    
    % Sampling
    DistanceSampling = SortNormDistance( 1 : limit - 1);
    score = zeros(1,5);
    breaks = [1 0.8 0.6 0.4 0.2 0];
    
    for i = 1 : 5
        
        dist1 = DistanceSampling < breaks(i);
        dist2 = DistanceSampling > breaks(i + 1);
        dist = DistanceSampling(logical(dist1.*dist2));
        len = length(dist);
        score(i) = (len/(limit - 1)) ^ factor * mean(dist);
        
    end
    
    score = score / sum(score);
    score = [0, round(score * samples)];    
    len = zeros(6, 1);
    sscore = cumsum(score);
    sampled = zeros(1, sscore(end));
    
    for i = 1 : 5
        dist1 = DistanceSampling < breaks(i);
        dist2 = DistanceSampling > breaks(i + 1);
        dist = DistanceSampling(logical(dist1.*dist2));
        len(i + 1) = len(i) + length(dist);
        sampled(sscore(i) + 1 : sscore( i + 1 )) = randsample(len(i) + 1 : len( i + 1), score(i + 1), true);%muestrear con reemplazo   
    end
    
    SortNormDistance(:) = 0;
    SortNormDistance(sampled) = 1;
    NormDistance(idx) = SortNormDistance;
    Map2 = nan(size(ReadInfo.Map));
    Map2(~Indicator) = NormDistance;
    sampled = find(Map2(:) == 1);
    
    %Extract long and lat
    Lat = NaN(1, length(sampled));
    Long = Lat;
    for i = 1 : length(sampled)
        col = fix(sampled(i)/Dim1);  % integer part
        row = rem(sampled(i), Dim1);  % remainder
        
        if row > Dim1 || col > Dim2
            disp(i)
        end
        
        [Long(i), Lat(i)] = pix2map(R, row, col);%!!
    end
    
%     NewNanPositions = [isnan(Long)', isnan(Lat)'];
%     NewNanPositions = sum(NewNanPositions, 2);
%     NewNanPositions = logical(NewNanPositions);
%     
%     Long(NewNanPositions) = [];
%     Lat(NewNanPositions) = [];
            
    if show == 1
        clf
        geoshow(Map, R, 'DisplayType', 'surface');
        contourcmap('jet',0 : 0.05 : 1, 'colorbar', 'on', 'location', 'vertical')
        geoshow(Lat, Long, 'DisplayType', 'Point', 'Marker', 'o', ...
            'MarkerSize', 5, 'MarkerFaceColor', [.95 .9 .8], 'MarkerEdgeColor', ...
            'black', 'Zdata', 2 * ones(length(Long), 1));  
    end
    T = table(repmat(spName,length(Lat),1));
    %T = table('Size',[length(Lat), 1],'VariableTypes',{'string'});
    T.LAT = Lat';
    T.LONG = Long';
    T.Properties.VariableNames{1}='Name';
end
