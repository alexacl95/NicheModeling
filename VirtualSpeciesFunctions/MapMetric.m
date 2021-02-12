function Metric = MapMetric(NicheMap, ModelMap, show)
    
    % Deffining the Niche map without nans
    VectorNicheMap = NicheMap(:);
    IdxNicheMap = ~isnan(VectorNicheMap);
    VectorNicheMap = VectorNicheMap(IdxNicheMap);
    
    % Deffining the Model map without nans
    VectorModelMap = ModelMap(:);
    VectorModelMap = VectorModelMap(IdxNicheMap);
    
    % Sort the model map pixels according to the niche map values
    [VectorNicheMap, IdxModelMap] = sort(VectorNicheMap, 'ascend');
    VectorModelMap = VectorModelMap(IdxModelMap);
 
    VectorModelMap = cumsum(VectorModelMap, 'omitnan');
    VectorNicheMap = cumsum(VectorNicheMap, 'omitnan');
    
    % Maximum pixel value between both maps
    maxim = max(max(VectorNicheMap), max(VectorModelMap)); 
    
    % Normalization between 0-1 each sorted pixel
    VectorNicheMap = VectorNicheMap/maxim;
    VectorModelMap = VectorModelMap/maxim;
    
    DiffMap = abs(VectorNicheMap - VectorModelMap);
     
    LengDiffMap = length(DiffMap);
    LengDiffMap = linspace(0, 1, LengDiffMap);
    
    Metric = 1 - trapz(LengDiffMap, DiffMap);
    
    if show == 1         
        LengthNicheMap = length(VectorNicheMap);
        LengthNicheMap = linspace(0, 1, LengthNicheMap); 
        clf
        plot(LengthNicheMap, VectorModelMap, 'LineWidth', 2)
        hold on
        plot(LengthNicheMap, VectorNicheMap, 'LineWidth', 2)
        legend('Estimated','Original','Location','best')
        title(strcat(num2str(round(Metric*100,2)),'%'))
    end
    
end
