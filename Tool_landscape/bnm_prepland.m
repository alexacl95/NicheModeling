function out = bnm_prepland(Z,threshold,landsample,parallel)

if nargin<4
    parallel=false;
end
    
    indMap = Z(:,:,1);
    index = ~isnan(indMap(:));
    dim = size(Z);
    tsize = sum(index);
    N = round(tsize*landsample);
    samples = randsample(1:tsize, N, false);
    disp("--Sampling f(z) with "+num2str(N)+" samples--")
    T = zeros(N,dim(3));
    for i=1:dim(3)
        aux = Z(:,:,i);
        sampler = aux(index);
        T(:,i)=sampler(samples);
    end
    [Corr, ~] = corr(T, 'rows', 'complete');
    Corr = abs(Corr);
    
    Tdata = table();
    ex = [];
    ex2 = ex;
    vars = [];
    for i = 1:dim(3)
        if sum(i == ex) == 0
            if sum(Corr(i, :) > threshold) > 1
                f = find(Corr(i, :) > threshold);
                covM = cov(T{:, f});
                covMi = inv(covM);
                D1 = diag(covM);
                D2 = diag(covMi);
                bestReg = 1-1./(D1.*D2);
                [~,index] = sort(bestReg,'descend');
                regresor = f(index(1));
                ex = [ex, f];
                f = setdiff(f, regresor);
                ex2 = [ex2, f];
                %Kfold = 10;
                %LengthData = length(T{:, i});
                lambda_opt = k_fCV(T(:, i), T(:, f),20,parallel);
                %lambda_opt=5;
                %Revisar q pasa cuando no encuentra un lambda optimo
                if isempty(lambda_opt)
                    lambda_opt = 0.1;
                end                
                model = ridge([T{:,regresor}],[T{:,f}],lambda_opt,0);
                vars = [vars, regresor];
                Tdata = addprop(Tdata, {strcat('m', num2str(regresor))}, {'table'});
                eval(strcat("Tdata.Properties.CustomProperties.m", num2str(regresor), "=@(X) X(:,regresor)-(X(:,f)*model(2:end)+model(1));"))
            end
        end
    end
    
    disp('----Creating landscape predictors----')

    %indicators=setdiff(1:size(T,2),union(ex,vars));
    indicators = setdiff(1 : size(T, 2), ex);
    siz = size(indicators, 2);
    D1 = floor(sqrt(siz)); % Number of rows of subplot
    D2 = D1 + ceil((siz - D1^2)/D1);
    count = 1;
    
    for i = indicators
        [funcKernel, xdata] = ksdensity(T(:, i));    
        funcKernel = normalize(funcKernel, 'range');
        eval(strcat("Tdata.bio", num2str(i), "= [xdata; funcKernel]';"))
        count = count+1;
    end
    
    siz = size(vars,2);
    D1 = floor(sqrt(siz)); % Number of rows of subplot
    D2 = D1 + ceil((siz-D1^2)/D1);
    count = 1;

    for i = vars
        [funcKernel, xdata] = ksdensity(eval(strcat("Tdata.Properties.CustomProperties.m", num2str(i), "(T)")));
        funcKernel = normalize(funcKernel, 'range');
        eval(strcat("Tdata.var", num2str(i), "= [xdata; funcKernel]';"))
        count = count + 1; 
    end
    
    out = struct();
    out.Tdata = Tdata;
    out.Indicators = indicators;
    out.Vars = vars;
    
    disp('----Done!----')

end