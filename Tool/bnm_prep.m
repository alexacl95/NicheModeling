function outFinal = bnm_prep(doc,layerfolder,show,threshold,cluster)

tic

if nargin < 4
    threshold = 0.8;
end

if nargin < 5
    cluster = true;
end

if istable(doc)
    T = doc;
else
    T = readtable(doc);
end

if isstring(layerfolder) || ischar(layerfolder)
    temp = dir(layerfolder);
    layers = {temp.name};
    comp = strncmp('bio',layers,3);
    layers = layers(comp);
    disp('----Reading layers----')

    [Z,R] = arcgridread(strcat(layerfolder,layers{1}));
    
    N = length(layers);
    Z = zeros(size(Z,1),size(Z,2),N);

    for i = 1:N
%        progressbar(i/N)
       Z(:, :, i) = arcgridread(strcat(layerfolder,layers{i}));
       eval(strcat('T.bio',num2str(i),'=ltln2val(Z(:,:,i),R,[T.LAT],[T.LONG]);'))
    end
    toc
else
    
    Z = layerfolder.Z;
    R = layerfolder.R;
    N = size(Z,3);
    
    for i = 1:N
%        progressbar(i/N)
       eval(strcat('T.bio', num2str(i), '=ltln2val(Z(:,:,i),R,[T.LAT],[T.LONG]);'))
    end
end

if cluster
    %%% Finding clusters with spectral clustering %%%
    disp('----Finding clusters----')
    T = bnm_clustering(T);
    %toca arreglar el numero minimo de puntos por cluster
    %%% ----------------------------------------- %%%
    K = T.Properties.CustomProperties.NumClusters;
    ClusterIndex = T.Properties.CustomProperties.ClusterIndex;
    disp(num2str(K)+" clusters identified")
    
    if show
        figure
        colormap(bone)
        map = Z(:,:,1);
        map(map>0)=0;
        geoshow(map,R,'DisplayType','texturemap')
        hold on
        color = lines(K); % Generate color values
        gscatter(T.LONG,T.LAT, ClusterIndex,color(1:K,:));
        axis off
    end
else
    K=1;
    ClusterIndex = true(1,size(T,1));
end


outFinal = cell(1, length(K));


T_Original = T;

for ij = 1:K
    
    T = T_Original;
    indexClusterij = find(ClusterIndex == ij);
    
    disp('----Finding correlation----')
    

    %Construction of new variables
    
    T2 = T(indexClusterij,:);
    
    T = T(indexClusterij, 4:end);
    
    if show
        correlationCircles(T,'varnames',T.Properties.VariableNames(4:end))
    end
    
    [Corr, ~] = corr(T{:,:}, 'rows', 'complete');
    Corr = abs(Corr);
    
    Tdata = table();
    ex = [];
    vars = [];

    for i = 1:N
        if sum(i == ex) == 0
            if sum(Corr(i, :) > threshold) > 1
                f = find(Corr(i, :) > threshold);
                f = setdiff(f, i);
                ex = [ex, f];
                %Kfold = 10;
                %LengthData = length(T{:, i});
                lambda_opt = k_fCV([T{:, i}], [T{:, f}]);
                %Revisar q pasa cuando no encuentra un lambda optimo
                if isempty(lambda_opt)
                    lambda_opt = 0.1;
                end                
                model = ridge([T{:,i}],[T{:,f}],lambda_opt,0);
                vars = [vars, i];
                Tdata = addprop(Tdata, {strcat('m', num2str(i))}, {'table'});
                eval(strcat("Tdata.Properties.CustomProperties.m", num2str(i), "=@(X) X{:,i}-(X{:,f}*model(2:end)+model(1));"))
            end
        end
    end

    toc
    disp('----Creating predictors----')

    %indicators=setdiff(1:size(T,2),union(ex,vars));
    indicators = setdiff(1 : size(T, 2), ex);
    siz = size(indicators, 2);
    D1 = floor(sqrt(siz)); % Number of rows of subplot
    D2 = D1 + ceil((siz - D1^2)/D1);
    count = 1;

    %outlier remotion must starts before saving Tdata as predictors.
    %It is highly recommended to start outlier remotion from normalized PCA
    %table

    if show
        figure('Name', 'Histogram')
    end

    for i = indicators
        if show
            subplot(D1, D2, count)
            histfit(T{:, i}, 10, 'kernel');
            title(strcat('bio', num2str(i)))
        end
        [funcKernel, xdata] = ksdensity(T{:, i});    
        funcKernel = normalize(funcKernel, 'range');
        eval(strcat("Tdata.bio", num2str(i), "= [xdata; funcKernel]';"))
        count = count+1;
    end

    if show
        figure('Name','Histogram2')
    end

    siz = size(vars,2);
    D1 = floor(sqrt(siz)); % Number of rows of subplot
    D2 = D1 + ceil((siz-D1^2)/D1);
    count = 1;

    for i = vars
        if show
            subplot(D1,D2,count)
            h = histfit(eval(strcat("Tdata.Properties.CustomProperties.m", num2str(i), "(T)")), 10, 'kernel');
            title(strcat('var',num2str(i)))
        end
        [funcKernel, xdata] = ksdensity(eval(strcat("Tdata.Properties.CustomProperties.m", num2str(i), "(T)")));
        funcKernel = normalize(funcKernel, 'range');
        eval(strcat("Tdata.var", num2str(i), "= [xdata; funcKernel]';"))
        count = count + 1; 
    end
    
    out = struct();
    
    out.Tdata = Tdata;
    out.Indicators = indicators;
    out.Vars = vars;
    out.T2 = T2;
    out.Z = Z;
    out.R = R;
   
    outFinal{ij} = out;   
    
    
end

toc
    
disp('¡All done!')

end