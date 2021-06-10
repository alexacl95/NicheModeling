function T = samplingVS(ReadInfo, MapInfo, samples, factor, show, spName, replacement, uniqueness)

    if nargin<8
        uniqueness = false;
    end
    
    Map = MapInfo.Map;
    SortNormDistance = MapInfo.SortNormDistance;
    NormDistance = MapInfo.NormDistance;
    R = ReadInfo.R;
    idx = MapInfo.idx;
    Indicator = ReadInfo.Indicator;

% Sampling
    switch factor >= 0
        case true
            limit = find(SortNormDistance, 1, 'last')+1;
            DistanceSampling = SortNormDistance( 1 : limit-1);
            score = zeros(1,5);
            breaks = [1 0.8 0.6 0.4 0.2 0];

            if factor <=1
                for i = 1 : 5

                    dist1 = DistanceSampling < breaks(i);
                    dist2 = DistanceSampling > breaks(i + 1);
                    dist = DistanceSampling(logical(dist1.*dist2));
                    len = length(dist);
                    score(i) = (len/(limit - 1)) ^ factor * mean(dist);

                end
            else
                for i = 1 : 5

                    dist1 = DistanceSampling < breaks(i);
                    dist2 = DistanceSampling > breaks(i + 1);
                    dist = DistanceSampling(logical(dist1.*dist2));
                    len = length(dist);
                    if i<=2
                        score(i) = mean(dist);
                    else
                        score(i) = (len/(limit - 1)) * mean(dist);
                    end

                end
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
                sampled(sscore(i) + 1 : sscore( i + 1 )) =...
                    randsample(len(i) + 1 : len( i + 1), score(i + 1), replacement,...
                    DistanceSampling(len(i) + 1 : len( i + 1))/...
                    sum(DistanceSampling(len(i) + 1 : len( i + 1))));%muestrear   
            end
            
        case false
            
            DistanceSampling = 1:length(SortNormDistance);
            factor=abs(factor);
            if factor>1
                SortNormDistance=SortNormDistance.^factor;
            end
            sampled = randsample(DistanceSampling,round(samples),replacement,...
                SortNormDistance/sum(SortNormDistance));
            
            
        otherwise
            T = [];
            disp('Wrong factor, try "-1" or "0"')
            return
    end
    
    Map2 = nan(size(ReadInfo.Map));
    
    if ~uniqueness
   
        Tsamples = length(sampled);
        sampledRow = zeros(1,Tsamples);
        sampledCol = sampledRow;
        for i = 1:Tsamples
            SortNormDistance(:) = 0;
            SortNormDistance(sampled(i)) = 1;
            NormDistance(idx) = SortNormDistance;
            Map2(~Indicator) = NormDistance;
            [sampledRow(i) , sampledCol(i)] = find(Map2 == 1);
            Map2(~Indicator) = 0;
        end
        sampledRow = sampledRow';
        sampledCol = sampledCol';
    else
        
        SortNormDistance(:) = 0;
        SortNormDistance(sampled) = 1;
        NormDistance(idx) = SortNormDistance;
        Map2(~Indicator) = NormDistance;
        [sampledRow , sampledCol] = find(Map2 == 1);
    end
    
    
    %Extract long and lat
    
    [Long, Lat] = pix2map(R,sampledRow,sampledCol);
    
%     Lat = NaN(1, length(sampled));
%     Long = Lat;
%        
%     for i = 1 : length(sampled)
%         col = fix(sampled(i)/Dim1);  % integer part
%         row = rem(sampled(i), Dim1);  % remainder
%         
%         if row > Dim1 || col > Dim2
%             disp(i)
%         end
%         
%         [Long(i), Lat(i)] = pix2map(R, row, col);%!!
%     end
    
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
    T.LAT = Lat;
    T.LONG = Long;
    T.Properties.VariableNames{1}='Name';
end