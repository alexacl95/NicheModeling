function T = samplingVS(ReadInfo, InfoInitialPoint, MapInfo, samples, factor, show, spName, replacement, uniqueness)
% T = samplingVS(ReadInfo, InfoInitialPoint, MapInfo, samples, factor, show,...
%                spName, replacement, uniqueness)
% 
% DESCRIPTION 
%   'samplingVS' sample the species presence-data(observations) over 
%   the generated map niche 
%
% REQUIRED INPUTS
%   ReadInfo: an strcuture generated by 'ReadLayers' function
%   InfoInitialPoint: and structure generated by 'InitialPoint' function
%   MapInfo: a structure generated by 'NicheGeneration' function
%   samples: a integer with the number of required samples 
%   factor (alpha): real number that define the sampling method
%   show: boolean variable (true, false), show the resulting niche map  
%         and the sampled points
%   spName: a string with the species name
%   replacement: boolean variable (true, false), perform a replacemet in 
%                the sample if true, or not replacement if false
%   
% OPTIONAL INPUTS
%   uniqueness: (default: false)
% 
% OUTPUTS:
%   T: a table with the species name 'spName', the longitude and 
%   the latitude of the virtual species observations
%%  
    if nargin < 9
        uniqueness = false;
    end
    
    Map = MapInfo.Map;
    SortNormDistance = MapInfo.SortNormDistance;
    NormDistance = MapInfo.NormDistance;
    R = ReadInfo.R;
    idx = InfoInitialPoint.idx;
    Indicator = ReadInfo.Indicator;

% Sampling the niche map
    switch factor >= 0
        %Positive alpha
        case true
            limit = find(SortNormDistance, 1, 'last') + 1;
            DistanceSampling = SortNormDistance( 1 : limit-1);
            score = zeros(1, 5);
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
                    if i <= 2
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
                
                sampled(sscore(i) + 1 : sscore( i + 1 )) = ...
                    randsample(len(i) + 1 : len( i + 1), score(i + 1), replacement,...
                    DistanceSampling(len(i) + 1 : len( i + 1))/...
                    sum(DistanceSampling(len(i) + 1 : len( i + 1))));%muestrear   
            end
        % Negative alpha    
        case false
            
            DistanceSampling = 1:length(SortNormDistance);
            factor = abs(factor);
            
            if factor > 1
                SortNormDistance = SortNormDistance.^factor;
            end 
            
            sampled = randsample(DistanceSampling, round(samples), replacement,...
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
    
    %Plot the niche map with the samples
    if show == 1
        clf
        geoshow(Map, R, 'DisplayType', 'surface');
        contourcmap('jet',0 : 0.05 : 1, 'colorbar', 'on', 'location', 'vertical')
        geoshow(Lat, Long, 'DisplayType', 'Point', 'Marker', 'o', ...
            'MarkerSize', 5, 'MarkerFaceColor', [.95 .9 .8], 'MarkerEdgeColor', ...
            'black', 'Zdata', 2 * ones(length(Long), 1));  
    end
    
    % OUTPUT STORAGE
    %Create the table with samples coordinates
    T = table(repmat(spName, length(Lat), 1));    
    T.LONG = Long;
    T.LAT = Lat;
    T.Properties.VariableNames{1} = 'Name';
end