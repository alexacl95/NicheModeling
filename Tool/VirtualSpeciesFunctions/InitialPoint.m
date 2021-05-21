function InfoInitialPoint = InitialPoint(ReadInfo, point, coeff)

    Dimension = ReadInfo.Dimensions(1);
    Dimension3 = ReadInfo.Dimensions(2);
    NormalizedClimVar = ReadInfo.NormalizedClimVar;
    Distance = zeros(1,Dimension);
    
    if nargin < 2
        point = rand(Dimension3, 1);
    end
    if nargin < 3
        coeff = rand(Dimension3, 1);
    end

    coeff = coeff/sum(coeff);
    point = point.*coeff;
    
    NormalizedClimVar = NormalizedClimVar.*coeff;     
    %point = datasample(NormalizedClimVar', 1)'; %Another point generation
    
    for i = 1 : Dimension
        Distance(i) = norm(point - NormalizedClimVar(:, i))...
                      * (2 - corr2(point, NormalizedClimVar(:, i)));
    end
    
    NormDistance = 1 - normalize(Distance, 2, 'range');
    [SortNormDistance, idx] = sort(NormDistance, 2, 'descend');
    
    InfoInitialPoint.NormDistance = NormDistance;
    InfoInitialPoint.idx = idx;
    InfoInitialPoint.SortNormDistance = SortNormDistance;
    InfoInitialPoint.coeff = coeff ;
    

end